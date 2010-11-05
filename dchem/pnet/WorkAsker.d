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

class WorkAskerGen:RemoteCCTask,Sampler {
    char[] connectionUrl;
    double maxDuration=525600.0;
    double maxDurationTot=525600.0;
    long maxEval=long.max;
    long maxWait=100;
    long ownerCacheSize=10;
    char[] precision="LowP";
    InputField evaluator;
    
    mixin(serializeSome("dchem.WorkAsker",`
    connectionUrl: the url to use to get the connection to the silos
    maxDuration: maximum duration in minutes for each context (one year: 525600)
    maxDurationTot: maximum duration in minutes for the workAsker (one year: 525600)
    maxEval: maximum number of evaluations
    maxWait: maximum number of retries waiting for a point to evaluate
    ownerCacheSize: amout of points whose owner is cached (10)
    precision: precision of the silos, either LowP or Real (LowP)
    evaluator: the mehod used to evaluate energy and forces`));
    mixin printOut!();
    mixin myFieldMixin!();
    WorkAsker!(LowP)[] wAskersLowP;
    WorkAsker!(Real)[] wAskersReal;
    
    bool verify(CharSink s){
        bool res=true;
        if (precision!="LowP" && precision!="Real"){
            log("precision should be either LowP or Real, not '")(precision)("' in field ")(myField)("\n");
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
        if (evaluator is null || (cast(Method)evaluator.contentObj)is null){
            dumper(s)("evaluator should be valid and of type Method in field")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    
    /// runs a work asker
    void run(LinearComm pEnv, CharSink log){
        auto m=cast(Method)evaluator.contentObj;
        assert(m!is null);
        m.setup(pEnv,log);
        auto maxTotTime=noToutWatcher.now()+maxDurationTot*60;
        while(true){
            auto cc=m.getCalculator(true,[]);
            if (cc is null) break;
            workOn(cc,maxTotTime);
        }
    }
    
    void workOn(CalculationContext cc){
        workOn(cc,noToutWatcher.now()+maxDurationTot*60);
    }
    /// starts the task with the given CalculationContext
    void workOn(CalculationContext cc,ev_tstamp maxTotTime){
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
        synchronized(this){
            size_t i=wAskersLowP.length;
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
                auto pos=find(w,wAskersLowP);
                if (pos!=wAskersLowP.length){
                    wAskersLowP[pos]=wAskersLowP[wAskersLowP.length-1];
                    wAskersLowP[wAskersLowP.length-1]=null;
                    wAskersLowP=wAskersLowP[0..wAskersLowP.length-1];
                }
            }
        }
        static if (is(T==Real)){
            synchronized(this){
                auto pos=find(w,wAskersReal);
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
    CalculationContext ctx;
    double timeEnd;
    long leftEvals;
    long maxWait;
    enum Status{
        Configure,Running,Stopping,Stopped
    }
    Status status;
    
    this(WorkAskerGen input,CalculationContext ctx){
        this.input=input;
        auto silosGen=new PNetSilosClientGen(connectionUrl,ownerCacheSize);
        this.silos=new PNetSilosClient(silosGen,ctx);
        this.ctx=ctx;
        timeEnd=input.maxDuration*60+ev_time();
        leftEvals=input.maxEval;
        maxWait=input.maxWait;
        stop=Status.Configure;
    }
    void run(){
        try{
            synchronized(this){
                if (status!=Status.Configure){
                    version(TrackWorkAsker) {
                        sinkTogether(sout,delegate void(CharSink s){
                            dumper(s)("failed start of ")(this)("\n");
                        });
                    }
                    if (status==Status.Running){
                        throw new Exception("run called on already running WorkAsker",__FILE__,__LINE__);
                    }
                    return;
                }
                status=Status.Running;
            }
            version(TrackWorkAsker) {
                sinkTogether(sout,delegate void(CharSink s){
                    dumper(s)("started ")(this)("\n");
                });
            }
            while(leftEvals>0 && timeEnd>ev_time() && status==Status.Running){
                auto p=silos.pointToEvaluate(SKeyVal.All);
                version(TrackWorkAsker) {
                    sinkTogether(sout,delegate void(CharSink s){
                        dumper(s)(this)(" should eval point ")(p)("\n");
                    });
                }
                if (p.data == 0){
                    break;
                }
                if (p.data == 1){
                    if (maxWait<0) break;
                    --maxWait;
                    continue;
                } else {
                    maxWait=input.maxWait;
                }
                auto time=ev_time();
                try{
                    auto localP=silos.createLocalPoint(p,time);
                    if (localP.evalWithContext(ctx)){
                        version(TrackWorkAsker) {
                            sinkTogether(sout,delegate void(CharSink s){
                                dumper(s)(this)(" did eval point ")(p)("\n");
                            });
                        }
                        --maxEval;
                    } else {
                        version(TrackWorkAsker) {
                            sinkTogether(sout,delegate void(CharSink s){
                                dumper(s)(this)(" did *not* eval point ")(p)("\n");
                            });
                        }
                    }
                } catch (Exception e){
                    silos.evaluationFailed(SKeyVal.All,p);
                }
            }
        } catch(Exception e){
            sinkTogether(sout,delegate void(CharSink s){
                dumper(s)("Exception in WorkAsker ")(this)(":")(e);
            });
        }
        version(TrackWorkAsker) {
            sinkTogether(sout,delegate void(CharSink s){
                dumper(s)("finished")(this)("\n");
            });
        }
        input.rmWorkAsker(this);
        status=Status.Stopped;
        ctx.giveBack();
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
