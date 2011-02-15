module dchem.pnet.WorkAsker;
import dchem.Common;
import dchem.calculator.CalculatorModels;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.parallel.mpi.MpiModels;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.parallel.smp.Smp;
import blip.serialization.Serialization;
import dchem.pnet.PNetSilosClient;
import blip.io.EventWatcher: ev_tstamp,ev_time;
import blip.io.Console;
import blip.parallel.rpc.Rpc;
import blip.stdc.stdlib: abort;
import blip.core.Array;

version(TrackPNet){
    version=TrackWorkAsker;
}

class WorkAskerGen:RemoteCCTask,Sampler,SilosConnectorI {
    string _connectionUrl;
    string mainWGenUrl;
    double maxDuration=525600.0;
    double maxDurationTot=525600.0;
    long maxEval=long.max;
    long maxWait=100;
    long ownerCacheSize=10;
    string _precision="LowP";
    InputField evaluator;
    bool _noWorkLeft=false;
    
    mixin(serializeSome("dchem.WorkAsker",`
    connectionUrl: the url to use to get the connection to the silos
    maxDuration: maximum duration in minutes for each context (one year: 525600)
    maxDurationTot: maximum duration in minutes for the workAsker (one year: 525600)
    maxEval: maximum number of evaluations
    maxWait: maximum number of retries waiting for a point to evaluate
    ownerCacheSize: amout of points whose owner is cached (10)
    precision: precision of the silos, either LowP or Real (LowP)
    evaluator: the mehod used to evaluate energy and forces (used only if you run this standalone)`));
    mixin printOut!();
    mixin myFieldMixin!();
    WorkAsker!(LowP)[] wAskersLowP;
    WorkAsker!(Real)[] wAskersReal;
    
    void noWorkLeft(){
        _noWorkLeft=true;
        if (mainVendor is null && mainWGenUrl.length>0){
            rpcManualVoidCall(mainWGenUrl~"/noWorkLeft");
        }
    }

    mixin(rpcMixin("dchem.WorkAsker", "","noWorkLeft"));
    DefaultVendor mainVendor;
    
    /// precision of the silos
    string precision(){
        return _precision;
    }
    void precision(string p){
        _precision=p;
    }
    /// connection url
    string connectionUrl(){
        return _connectionUrl;
    }
    /// sets the connectionUrl
    void connectionUrl(string c){
        _connectionUrl=c;
    }
    /// sets precision and url
    void setConnectionAndPrecision(string c,string p){
        _precision=p;
        _connectionUrl=c;
    }
    
    bool verify(CharSink s){
        bool res=true;
        auto log=dumper(s);
        /+if (connectionUrl.length==0){
            log("connectionUrl has to be given to WorkAsker ")(myFieldName)("\n");
            res=false;
        }+/
        if (precision!="LowP" && precision!="Real"){
            log("precision should be either LowP or Real, not '")(precision)("' in field ")(myFieldName)("\n");
            res=false;
        }
        if (maxDuration<=0){
            dumper(s)("maxDuration should be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        if (maxDurationTot<=0){
            dumper(s)("maxDurationTot should be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        if (maxEval<=0){
            dumper(s)("maxEval should be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        if (evaluator !is null && (cast(Method)evaluator.contentObj)is null){
            dumper(s)("evaluator if given should be of type Method in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    
    /// runs a work asker
    void run(LinearComm pEnv, CharSink log){
        assert(pEnv!is null,"pEnv");
        assert(log!is null,"log");
        if (connectionUrl.length==0){
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("connectionUrl has to be given to WorkAsker ")(myFieldName)("\n");
            });
        }
        if (mainWGenUrl.length==0){
            if (mainVendor is null){
                if (pEnv.myRank==0){
                    mainVendor=new DefaultVendor(this);
                    ProtocolHandler.defaultProtocol.publisher.publishObject(mainVendor,"WorkAsker_"~myFieldName,true);
                    mainWGenUrl=mainVendor.proxyObjUrl();
                }
                mpiBcastT(pEnv,mainWGenUrl,0);
            }
        }
        auto m=cast(Method)evaluator.contentObj;
        assert(m!is null);
        m.setup(pEnv,log);
        auto maxTotTime=ev_time()+maxDurationTot*60;
        while(!_noWorkLeft){
            auto cc=m.getCalculator(true,[]);
            if (cc is null) break;
            cc.executeLocally(this);
        }
    }
    
    void workOn(LocalCalculationContext cc){
        workOn(cc,ev_time()+maxDurationTot*60);
    }
    /// starts the task with the given CalculationContext
    void workOn(LocalCalculationContext cc,ev_tstamp maxTotTime){
        if (connectionUrl.length==0){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("connectionUrl has to be given to WorkAsker ")(myFieldName)("\n");
            }),__FILE__,__LINE__);
        }
        switch(precision){
        case "LowP":
            auto lp=new WorkAsker!(LowP)(this,cc);
            if (lp.timeEnd>maxTotTime) lp.timeEnd=maxTotTime;
            synchronized(this){
                wAskersLowP~=lp;
            }
            Task("WorkOnStart",&lp.start).autorelease.submitYield();
            break;
        case "Real":
            auto rp=new WorkAsker!(Real)(this,cc);
            if (rp.timeEnd>maxTotTime) rp.timeEnd=maxTotTime;
            synchronized(this){
                wAskersReal~=rp;
            }
            Task("WorkOnStart",&rp.start).autorelease.submitYield();
            break;
        default:
            assert(0,"invalid precision "~precision);
        }
    }
    /// might stop the task, or not, might return immediatly (even if the task is still running)
    void stop(){
        size_t i;
        synchronized(this){
            i=wAskersLowP.length;
        }
        while(i!=0){
            WorkAsker!(LowP) w;
            synchronized(this){
                --i;
                if (wAskersLowP.length<=i) continue;
                w=wAskersLowP[i];
            }
            w.stop();
        }
        synchronized(this){
            i=wAskersReal.length;
        }
        while(i!=0){
            WorkAsker!(Real) w;
            synchronized(this){
                --i;
                if (wAskersReal.length<=i) continue;
                w=wAskersReal[i];
            }
            w.stop();
        }
    }
    void rmWorkAsker(T)(WorkAsker!(T) w){
        static if (is(T==LowP)){
            synchronized(this){
                auto pos=find(wAskersLowP,w);
                if (pos!=wAskersLowP.length){
                    wAskersLowP[pos]=wAskersLowP[wAskersLowP.length-1];
                    wAskersLowP[wAskersLowP.length-1]=null;
                    wAskersLowP=wAskersLowP[0..wAskersLowP.length-1];
                }
            }
        }
        static if (is(T==Real)){
            synchronized(this){
                auto pos=find(wAskersReal,w);
                if (pos!=wAskersReal.length){
                    wAskersReal[pos]=wAskersReal[wAskersReal.length-1];
                    wAskersReal[wAskersReal.length-1]=null;
                    wAskersReal=wAskersReal[0..wAskersReal.length-1];
                }
            }
        }
    }
}

class WorkAsker(T){
    WorkAskerGen input;
    PNetSilosClient!(T) silos;
    LocalCalculationContext ctx;
    double timeEnd;
    long leftEvals;
    long maxWait;
    enum Status{
        Configure,Running,Stopping,Stopped
    }
    Status status;
    
    void logMsg(void delegate(CharSink)writer){
        silos.logMsg(delegate void(CharSink s){
            dumper(s)("<WorkAskerLog T=")(T.stringof)(" field=\"")(input.myFieldName)("\" addr=")(cast(void*)this)(">\n  ");
            indentWriter("  ",s,writer);
            s("\n</WorkAskerLog>");
        });
    }
    
    this(WorkAskerGen input,LocalCalculationContext ctx){
        this.input=input;
        auto silosGen=new PNetSilosClientGen(input.connectionUrl,input.ownerCacheSize);
        this.silos=new PNetSilosClient!(T)(silosGen,ctx);
        this.ctx=ctx;
        timeEnd=input.maxDuration*60+ev_time();
        leftEvals=input.maxEval;
        maxWait=input.maxWait;
        status=Status.Configure;
    }
    void start(){
        logMsg(delegate void(CharSink s){ s("starting work asker"); });
        assert(status<=Status.Running);
        bool noWorkLeft=false;
        scope(exit){
            input.rmWorkAsker(this);
            status=Status.Stopped;
            ctx.giveBack();
        }
        try{
            synchronized(this){
                if (status!=Status.Configure){
                    logMsg(delegate void(CharSink s){ dumper(s)("failed start"); });
                    if (status==Status.Running){
                        throw new Exception("run called on already running WorkAsker",__FILE__,__LINE__);
                    }
                    return;
                }
                status=Status.Running;
            }
            version(TrackWorkAsker) {
                logMsg(delegate void(CharSink s){
                    dumper(s)("started ")(this);
                });
            }
            while(leftEvals>0 && timeEnd>ev_time() && status==Status.Running){
                auto newOp=silos.getNextOp(SKeyVal.Master);
                if (newOp!is null){
                    newOp.initOp(ctx,silos);
                    version(TrackWorkAsker) {
                        logMsg(delegate void(CharSink s){
                            dumper(s)("should eval ")(newOp);
                        });
                    }
                    --leftEvals;
                    newOp.doOp();
                } else {
                    input.noWorkLeft();
                    break;
                }
            }
        } catch(Exception e){
            logMsg(delegate void(CharSink s){
                dumper(s)("Exception in WorkAsker:")(e);
            });
            abort();
        }
        version(TrackWorkAsker) {
            logMsg(delegate void(CharSink s){
                dumper(s)("finished");
            });
        }
    }
    
    void stop(){
        synchronized(this){
            if (status>=Status.Stopping){
                return;
            }
            status=Status.Stopping;
        }
        ctx.stop();
    }
}
