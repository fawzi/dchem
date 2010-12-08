module dchem.pnet.PNetSilos;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import dchem.input.RootInput;
import dchem.calculator.Calculator;
import dchem.calculator.CalculatorModels;
import dchem.pnet.PNetModels;
import dchem.pnet.MainPoint;
import blip.core.Variant;
import blip.container.HashMap;
import blip.sync.UniqueNumber;
import blip.sync.Atomic;
import blip.container.MinHeap;
import blip.container.GrowableArray;
import blip.container.Set;
import blip.math.random.Random;
import dchem.pnet.DenseLocalPointArray;
import dchem.sys.Constraints;
import blip.util.NotificationCenter;
import dchem.input.WriteOut;
import blip.math.IEEE;
import blip.time.Clock;
import blip.time.Time;
import blip.parallel.smp.WorkManager;
import blip.container.Deque;
import blip.container.Cache;
import blip.container.Pool;
import dchem.sys.DynVars;
import blip.util.LocalMem;
import blip.parallel.smp.Wait;
import blip.io.Console;
import blip.parallel.rpc.RpcMixins;
import blip.parallel.mpi.Mpi;
import blip.io.EventWatcher;
//import dchem.pnet.WorkAsker;

/// help structure 
struct CachedPoint(T){
    MainPoint!(T) mainPoint;
    ev_tstamp lastSync;
}

struct LoadStats{
    Real load;
    SKey silosKey;
    mixin(serializeSome("LoadStats","load|silosKey"));
    mixin printOut!();
    int opCmp(LoadStats l2){
        return ((load<l2.load)?-1:((load==l2.load)?cmp(cast(long)silosKey,cast(long)l2.silosKey):1));
    }
}

class SilosRegistry(T){
    HashMap!(char[],LocalSilosI!(T)) localSilos;
    void registerSilos(char[]url,LocalSilosI!(T) silos){
        synchronized(localSilos){
            auto lSil=url in localSilos;
            if (lSil !is null && silos!is (*lSil)){
                throw new Exception("overwriting silos for url "~url,__FILE__,__LINE__);
            }
            localSilos[url]=silos;
        }
    }
    LocalSilosI!(T) proxyForUrl(char[] url){
        synchronized(localSilos){
            auto lSil=url in localSilos;
            if (lSil !is null){
                (*lSil).retain;
                return lSil;
            }
        }
        
        auto prx=ProtocolHandler.proxyForUrl(url);
        auto silos=cast(PNetSilosI!(T))cast(Object)prx;
    }
}


class RemotePointEval(T):RemoteCCTask{
    Point point;
    char[] silosCoreUrl;
    CalculationContext ctx;
    LocalSilosI!(T) localSilos;
    mixin(serializeSome("dchem.minEE.RemotePointEval!("~T.stringof~")","point|silosCoreUrl"));
    
    this(){}
    this(Point p,char[]oUrl){
        point=p;
        silosCoreUrl=oUrl;
    }
    void workOn(CalculationContext ctx){
        //localSilos=
        //auto ctx=cast(CalculationContext)cast(Object)ProtocolHandler.proxyForUrl(ownerUrl);
        //auto pPos=ctx.pointPos(point); // get real point
        //ctx=args.get!(CalculationContext)();
    }
    void stop(){
        if (ctx!is null) ctx.stop();
        ctx=null;
    }
    
}

class SilosGen:Sampler{
    PNetSilos!(Real)[] silosReal;
    PNetSilos!(LowP)[] silosLowP;
    
    /// number of exploration steps to perform
    long explorationSteps;
    /// number of exploration steps to perform
    long maxExplorationTries=10;
    /// the current method for energy evaluation
    InputField evaluator;
    /// the task to execute on the evaluator (normally work asker)
    InputField evaluatorTask;

    /// base name of the silos
    char[] baseName="Silos";
    /// precision with which the points are stored
    char[] precision="LowP";
    /// list of observers and explorers
    InputField[] explorers;
    /// list of loaders (called after the explorer are set up, before exploring)
    InputField[] loaders; // use simply SilosWorkerGen[]??
    /// list of monitors, called at regular intervals during the exploration
    InputField[] monitors; // use simply SilosWorkerGen[]??
    /// list of finishers, called after the exploration is finished
    InputField[] finishers; // use simply SilosWorkerGen[]??
    /// the operations to execute on the contexts 

    // topology parameters
    /// minimum residual after projecting out the linear dependend directions
    Real minProjectionResidual=0.05;
    /// minimum cos of then angle in dual space to consider a direction equivalent
    Real sameDirCosAngle=0.866;
    /// minimum norm in the dual T space to accept an exploration (in units of explorationSize)
    Real minNormDual=0.4;
    /// minimum norm in the dual T space to accept an exploration for a self generated direction 
    /// (in units of discretizationStep)
    Real minNormDualSelf=0.1;
    /// discretization step
    Real discretizationStep=0.2;
    /// minimum cartesian difference
    Real cartesianDiffMin=0.02;
    /// max norm in the dual T space to accept an exploration (in units of explorationSize) defines the neighbors
    Real maxNormDual=1.9;
    /// square of the maximum distance from the optimal direction point in explorationSize units
    Real dirSize2=1.0;
    /// distance that is considered equivalent to 0
    Real zeroLen=1.e-10;
    
    // implicit topology parameters
    /// minimum norm in the real (cartesian) space for a self generated direction to be accepted before moving
    Real minRealNormSelf0(){
        return cartesianDiffMin/4;
    }
    /// minimum norm in the real (cartesian) space to which a self generated direction should be rescaled
    Real minRealNormSelf1(){
        return cartesianDiffMin;
    }
    /// minimum norm in the real (cartesian) space for the real movement (after constraints,...)
    /// of a self generated direction to be accepted
    Real minRealNormSelf2(){
        return cartesianDiffMin/2;
    }
    /// explorationSize
    Real explorationSize(){
        return discretizationStep;
    }
    /// square of maximum distance from the optimal direction point in cartesian units
    Real dirCartesianSize2(){
        return cartesianDiffMin*cartesianDiffMin;
    }
    Real maxMoveInternal(){
        return 0.0;
    }
    Real maxMoveCartesian=0.0;
    Real maxMoveDual=3.0;
    Real minMoveInternal=0.0;
    Real maxNormCartesian(){
        return 2.0*cartesianDiffMin;
    }
    Real sameDirCosAngleCartesian=0.866;
    Real dirDualSize2(){
        return 1.0;
    }
    Real inDirCartesianScale2=0.5;
    
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.Silos",`
        evaluator: the method to perform energy evaluations
        evaluatorTask: the task to execute on the evaluator (with no task you should start them manually)
        precision: the precision with which the points are stored
        explorers: the observers and explorers active on the point network
        loaders: SilosWorkers that load already evaluated points (called after the explorer are set up, before exploring)
        monitors: SilosWorkers called at regular intervals during the exploration
        finishers: list of finishers, SilosWorkers called after the exploration is finished
        minProjectionResidual: minimum residual after projecting out the linear dependend directions (0.05)
        sameDirCosAngle: minimum cos of then angle in dual space to consider a direction equivalent (0.866)
        minNormDual: minimum norm in the dual T space to accept an exploration (in units of explorationSize) (0.4)
        minNormDualSelf: minimum norm in the dual T space to accept an exploration for a self generated direction (in units of discretizationStep) (0.1)
        discretizationStep: discretization step (0.2)
        cartesianDiffMin: minimum cartesian difference (0.02)
        maxNormDual: max norm in the dual T space to accept an exploration (in units of explorationSize) (1.9)
        dirSize2: square of the maximum distance from the optimal direction point in explorationSize units (1.0)
        explorationSteps: number of exploration steps to perform (0)
        maxExplorationTries: number of times an exploration attempts to get a vaid point to evaluate before giving up (10)
        inDirCartesianScale2: square of the scale applied for the error in cartesian units on the component in the direction (0.5)
        sameDirCosAngleCartesian: cos of the angle in cartesian units (0.866)
        maxMoveCartesian: quick cutoff for point further than this in cartesian units (0.0)
        maxMoveDual: quick cutoff for points further than this in cartesian units (3.0)
        minMoveInternal: minimum difference in the internal coords to accept a new direction (0.0)
        zeroLen: distance that is considered equivalent to 0 (1.0e-10)
    `));
    this(){}

    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        if (evaluator is null || cast(Method)evaluator.contentObj is null){
            log("evaluator should be valid and be a Method in field ")(myFieldName)("\n");
            res=false;
        }
        if (evaluatorTask !is null && cast(RemoteCCTask)evaluatorTask.contentObj is null){
            log("evaluatorTask if given should be a RemoteCCTask in field ")(myFieldName)("\n");
            res=false;
        }
        if (precision!="LowP" && precision!="Real"){
            log("precision should be either LowP or Real, not '")(precision)("' in field ")(myFieldName)("\n");
            res=false;
        }
        foreach(e;explorers){
            auto o=cast(Object)e.content;
            if (cast(ExplorationObserverGen)o is null){
                log("explorers should be of either explorers or observer, not '")(o.classinfo.toString())("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(l;loaders){
            auto o=cast(Object)l.content;
            if (cast(SilosWorkerGen)o is null){
                log("loaders should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(m;monitors){
            auto o=cast(Object)m.content;
            if (cast(SilosWorkerGen)o is null){
                log("monitors should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(f;finishers){
            auto o=cast(Object)f.content;
            if (cast(SilosWorkerGen)o is null){
                log("finishers should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        return res;
    }

    void run(LinearComm pWorld,CharSink log){
        ulong n;
        n=rand.uniformR2(cast(ulong)SKeyVal.FirstValid,ulong.max); // hopefully unique
        char[] name=collectAppender(delegate void(CharSink s){
            dumper(s)(baseName)(n);
        });
        SKey newKey=cast(SKey)n;
        if (precision=="LowP"){
            auto silos=new PNetSilos!(LowP)(this,name,newKey);
            synchronized(this){
                silosLowP~=silos;
            }
            silos.run(pWorld,log);
        } else if (precision=="Real"){
            auto silos=new PNetSilos!(Real)(this,name,newKey);
            synchronized(this){
                silosReal~=silos;
            }
            silos.run(pWorld,log);
        } else {
            assert(0,"invalid precision value");
        }
    }
    
    void stop(){
        synchronized(this){
            foreach (s;silosLowP){
                s.stop();
            }
            foreach (s;silosReal){
                s.stop();
            }
        }
    }
    
}

char[] realFromInput(char[]props){
    char[] res;
    auto propsP=ctfeSplit("| \n",props,true);
    foreach(property;propsP){
        if (property.length>0){
            res~=`
    T `~property~`(){
        return cast(T)input.`~property~`;
    }`;
        }
    }
    res~=`
    Real[char[]] propertiesDict0(){
        Real[char[]] res;`;
    foreach(property;propsP){
        if (property.length>0){
            res~=`
        res["`~property~`"]=input.`~property~`;`;
        }
    }
    res~=`return res;
    }`;
    return res;
}

/// helper class to track publishing of points
/// we hack in the reuse because we know that there is *always* one (and only one) waiter
class PointBcastProgress(T){
    PNetSilos!(T) localContext;
    MainPointI!(T) pointToBCast;
    HashMap!(SKey,int) publishLevel;
    long waitingName=1;
    WaitConditionT!(false) nameKnown;
    PointBcastProgress nextPoint;
    PoolI!(PointBcastProgress) pool;
    this(PoolI!(PointBcastProgress)p){ pool=p; }
    this(MainPointI!(T) p,PoolI!(PointBcastProgress)pool=null){
        pointToBCast=p;
        this.pool=pool;
        this.localContext=cast(PNetSilos!(T))cast(Object)p.localContext;
        this.nameKnown=new WaitConditionT!(false)(&this.nameIsKnown);
        assert(this.localContext!is null);
    }
    void clear(){
        pointToBCast=null;
    }
    static PoolI!(PointBcastProgress) gPool;
    static this(){
        gPool=cachedPool(function PointBcastProgress(PoolI!(PointBcastProgress)p){
            return new PointBcastProgress(p);
        });
    }
    void reset(MainPointI!(T) p){
        pointToBCast=p;
        this.localContext=cast(PNetSilos!(T))cast(Object)p.localContext;
        assert(this.localContext!is null);
    }
    static PointBcastProgress opCall(MainPointI!(T)p){
        auto nP=gPool.getObj();
        nP.reset(p);
        return nP;
    }
    void publish(){
        foreach(sK,sL;localContext.silos){
            synchronized(this){
                ++waitingName;
            }
            sL.publishPoint(sK,localContext._key,pointToBCast.point,pSysWriter(pointToBCast.pos),
                pointToBCast.explorationSize,pointToBCast.gFlags);
        }
        synchronized(this){
            --waitingName;
        }
    }
    bool nameIsKnown(){
        assert(waitingName>=0);
        return waitingName==0;
    }
    void communicateCollision(){
        bool addNew=false;
        auto dP=pointToBCast.driedPoint();
        synchronized(this){
            if (pointToBCast.drop){
                addNew=true;
            }
        }
        PointBcastProgress nextPP;
        if (addNew){
            auto np=localContext.newPointAt(pointToBCast.pos.dynVars.x);
            np[]=dP;
            nextPP=PointBcastProgress(np);
            synchronized(localContext.publishingPoints){
                localContext.publishingPoints[np.point]=nextPP;
            }
            nextPP.publish();
        }
        synchronized(this){
            if (nextPP!is null) nextPoint=nextPP;
            --waitingName;
        }
        nameKnown.checkCondition();
    }
    void communicatePublish(SKey origin){
        synchronized(this){
            --waitingName;
        }
        nameKnown.checkCondition();
    }
    MainPointI!(T) finalName(){
        nameKnown.wait();
        PointBcastProgress nextP;
        synchronized(this){
            nextP=nextPoint;
        }
        auto res=pointToBCast;
        if (nextP!is null) {
            res=nextP.finalName;
        }
        if (pool!is null) pool.giveBack(this);
        synchronized(localContext.publishingPoints){
            localContext.publishingPoints.removeKey(pointToBCast.point);
        }
        return res;
    }
}

class PendingEvals(T) {
    PNetSilos!(T) silos;
    HashMap!(char[],EvalOp!(T)) pendingOps;
    long nPending;
    WaitCondition waitPending;
    UniqueNumber!(ulong) nextLocalId;
    
    this(PNetSilos!(T) silos){
        this.silos=silos;
        this.pendingOps=new HashMap!(char[],EvalOp!(T))();
        this.waitPending=new WaitCondition(&this.hasNoPendingOps);
        this.nextLocalId=UniqueNumber!(ulong)(1);
    }
    bool hasNoPendingOps(){
        return nPending==0;
    }
    
    void addPendingOp(EvalOp!(T) op){
        auto id=collectAppender(delegate void(CharSink s){
            dumper(s)(silos._key)("-")(nextLocalId.next());
        });
        op.initMaster(id,silos,silos._key);
        atomicAdd(nPending,1);
        synchronized(this){
            pendingOps[op.id]=op;
        }
    }
    
    void updateEvalStatus(char[] opId,ModifyEvalOp!(T) op,EvalOp!(T).Status newStat){
        EvalOp!(T) opAtt;
        synchronized(pendingOps){
            auto pOp=opId in pendingOps;
            if (pOp!is null) {
                opAtt=*pOp;
            }
        }
        bool resubmitted=false;
        if (opAtt!is null){
            if (op !is null){
                op.modifyEvalOp(opAtt);
            }
            if (opAtt.updateStatus(newStat) && newStat==EvalOp!(T).Status.Failure){
                silos.addEvalOp(SKeyVal.Master,opAtt);
                resubmitted=true;
            }
        }
        switch (newStat){
            case EvalOp!(T).Status.ToDo:
                assert(0);
            case EvalOp!(T).Status.InProgress:
                break;
            case EvalOp!(T).Status.Failure:
            case EvalOp!(T).Status.Success:
                if (opAtt && !resubmitted){
                    synchronized(pendingOps){
                        pendingOps.removeKey(opId);
                    }
                    atomicAdd(nPending,-1);
                    if (silos.paraEnv.myRank!=0){
                        silos.updateEvalStatus(SKeyVal.Master,opId,null,newStat);
                    }
                    silos.waitPendingEnd.checkCondition();
                } else if (silos.paraEnv.myRank==0){
                    atomicAdd(nPending,-1);
                    silos.waitPendingEnd.checkCondition();
                }
                break;
            default:
                assert(0);
        }
    }
}

class PNetSilos(T): LocalSilosI!(T){
    RunLevel runLevel;
    SilosGen input;
    UniqueNumber!(ulong) nextLocalId;
    UniqueNumber!(ulong) nextPntNr;
    UniqueNumber!(uint)  nextTag;
    string nameId;
    SKey _key;
    bool hasKey(SKey k){ return _key==k||k==SKeyVal.Any||(k==SKeyVal.Master && paraEnv.myRank==0); }
    HashMap!(SKey,PNetSilosI!(T)) silos; // contains self
    MinHeapSync!(LoadStats) loads; 
    BatchedGrowableArray!(Point,batchSize) localPointsKeys; // indexes for local owned points
    HashMap!(Point,SKey) owner;
    HashMap!(Point,CachedPoint!(T)) localCache;
    HashMap!(Point,PointBcastProgress!(T)) publishingPoints;
    PendingEvals!(T) pendingEvals;
    DenseLocalPointArray!(MainPointI!(T)) localPoints; /// local points (position,...)
    string name(){
        return nameId;
    }
    /// random number generator
    RandomSync rand(){ return _rand; }
    RandomSync _rand;
    /// reference position/particle system
    ParticleSys!(T) refPos() { return _refPos; }
    ParticleSys!(T) _refPos;
    /// constraints (taken from the evaluator, guaranteed to be non null)
    ConstraintI!(T) _constraints;
    /// place to log messages
    CharSink log;
    /// notification center
    NotificationCenter nCenter(){
        return _nCenter;
    }
    NotificationCenter _nCenter;
    /// explorers
    HashMap!(char[],ExplorerI!(T)) explorers;
    Deque!(char[]) activeExplorers;
    WaitCondition waitExplorers;
    WaitCondition waitPendingEnd;
    Deque!(EvalOp!(T)) localOps;
    long nPendingOps;
    /// parallel enviroment of silos
    LinearComm silosParaEnv;
    /// observers
    Deque!(ExplorationObserverI!(T)) observers;
    /// evaluator
    Method evaluator;
    RemoteCCTask evaluatorTask;
    PoolI!(MainPoint!(T)) pPool;
    Deque!(SilosWorkerI!(T)) loaders;
    Deque!(SilosWorkerI!(T)) monitors;
    Deque!(SilosWorkerI!(T)) finishers;
    long nEvals;
    SKey[] sKeys; /// keys of the various core silos (index is their rank)
    
    MainPoint!(T)allocPoint(PoolI!(MainPoint!(T))p){
        return new MainPoint!(T)(this,p);
    }
    bool hasActiveExpl(){
        return activeExplorers.length>0 || explorers.length==0;
    }
    bool isStopping(){
        return runLevel>RunLevel.WaitPending;
    }
    /// writes out a log message
    void logMsg(void delegate(CharSink)writer){
        char[512] buf;
        auto gArr=lGrowableArray(buf,0);
        auto s=dumper(&gArr.appendArr);
        s("<MinEELog id=\"")(_key)("\" time=\"")(Clock.now.ticks)("\">");
        writer(s.call);
        s("</MinEELog>\n");
    }
    /// writes out a log message
    void logMsg(char[]msg){
        logMsg(delegate void(CharSink s){ s(msg); });
    }
    /// constraints of the current system
    ConstraintI!(T) constraints(){
        return _constraints;
    }
    /// if the gradient is cheap to compute
    bool cheapGrad(){
        return false;
    }
    /// if the gradient should be speculatively calculated. calls addEnergyEvalLocal (in background)
    bool speculativeGradient(SKey s,Point p,Real energy){
        bool res=false;
        int i=0;
        synchronized(explorers){
            foreach(k,e;explorers){
                res=res || e.speculativeGradientLocal(p,energy);
            }
        }
        struct EAdd{
            void delegate(SKey,Point,Real) op;
            SKey s;
            Point p;
            Real energy;
            void doOp(){
                op(s,p,energy);
            }
        }
        auto eAdd=new EAdd;
        eAdd.s=s;
        eAdd.energy=energy;
        eAdd.p=p;
        eAdd.op=&addEnergyEvalLocal;
        Task("addOp",&eAdd.doOp).autorelease.submit(defaultTask);
        return res;
    }
    
    mixin(realFromInput(propertiesList));
    
    Real[char[]] propertiesDict(SKey s){
        if (hasKey(s)){
            auto res=propertiesDict0();
            res["cheapGrad"]=cheapGrad();
            return res;
        } else {
            return silosForKey(s).propertiesDict(s);
        }
    }
    
    /// returns a globally unique string 
    char[] nextUniqueStr(){
        return collectAppender(delegate void(CharSink s){
            s(nameId); s("_"); writeOut(s,nextLocalId.next());
        });
    }
    
    bool isAtPendingEnd(){
        bool res=false,checkPending=false;
        synchronized(pendingEvals){
            volatile bool hasOps=pendingEvals.hasNoPendingOps;
            if ((hasOps && runLevel>RunLevel.Running) || runLevel>RunLevel.WaitPending){
                if (runLevel==RunLevel.WaitPending && paraEnv.myRank==0) {
                    checkPending=true;
                    runLevel=RunLevel.Stopping;
                }
                res=true;
            }
        }
        if (checkPending){
            waitPendingEnd.checkCondition();
        }
        return res;
    }
    mixin(descSome("dchem.PNetSilos_"~T.stringof,`_key|runLevel|nextLocalId|nextPntNr|nameId`));
    mixin printOut!();
    
    void startAskers(){
        if (evaluatorTask is null) return;
        try{
            /// start the askers
            while(true){
                CalculationContext cInstance=evaluator.getCalculator(true,[]);
                if (cInstance is null) break;
                cInstance.executeLocally(evaluatorTask);
            }
        } catch (Exception e){
            sinkTogether(log,delegate void(CharSink s){
                dumper(s)("error while starting askers:")(e)("\n");
            });
        }
    }
    
    void start(LinearComm pEnv, CharSink log){
        pEnv.barrier();
        silosParaEnv=pEnv;
        evaluator.setup(pEnv,log);
        CalculationContext cInstance;
        if (_refPos is null || _constraints is null){
            cInstance=evaluator.getCalculator(true,[]);
        }
        // setup refPos
        if (_refPos is null){
            switch(cInstance.activePrecision()){
            case Precision.Real:
                _refPos=cInstance.refPSysReal.dupT!(T)(PSDupLevel.All);
                break;
            case Precision.LowP:
                _refPos=cInstance.refPSysLowP.dupT!(T)(PSDupLevel.All);
                break;
            default:
                assert(0,"unknown activePrecision");
            }
        }
        if (_constraints is null){
            auto cGen=cInstance.constraintGen;
            if (cGen is null) {
                _constraints=new NoConstraint!(T)();
            } else {
                _constraints=constraintT!(T)(cGen,_refPos);
            }
        }
        // run the loaders
        foreach(sw;loaders){
            sw.workOn(this);
        }
        runLevel=RunLevel.Running;
        if (explorers.length==0){
            log("no explorers, stopping (to avoid this add a dchem.WaitExplorer)\n");
            increaseRunLevel(SKeyVal.All,RunLevel.WaitPending);
        }
        pEnv.barrier();
        sKeys=new SKey[](pEnv.dim);
        mpiAllGatherT(pEnv,_key,sKeys);
        if (sKeys.length>1){
            auto sortedK=sKeys.dup;
            sort(sortedK);
            foreach(i,v;sortedK[1..$]){
                if (sortedK[i]==v){
                    throw new Exception("key collision",__FILE__,__LINE__);
                }
            }
        }
        /// start the askers
        Task("startAskers",&startAskers).autorelease.submitYield();
    }
    void run(LinearComm pEnv, CharSink log){
        start(pEnv,log);
        waitPendingEnd.wait();
        if (paraEnv.myRank==0){
            increaseRunLevel(SKeyVal.All,RunLevel.Stopping);
        }
        pEnv.barrier();
        synchronized(pendingEvals){
            if (runLevel>RunLevel.Stopping) return;
            runLevel=RunLevel.Stopping;
        }
        // run the finishers
        foreach(sw;finishers){
            sw.workOn(this);
        }
        pEnv.barrier();
        runLevel=RunLevel.Stopped;
    }
    
    void activateExplorer(SKey key,char[] name)
    in{
        synchronized(explorers){
            assert(name in explorers);
        }
    } body {
        if (hasKey(key)){
            if (_key!=sKeys[0]) throw new Exception("unexpected sKey in activateExplorer",__FILE__,__LINE__);
            foreach(expl;activeExplorers){
                if (expl == name) return;
            }
            if (paraEnv.myRank==0){
                activeExplorers.pushBack(name);
                waitExplorers.checkCondition();
            } else {
                activateExplorer(SKeyVal.Master,name);
            }
        } else {
            silosForKey(key).activateExplorer(key,name);
        }
    }
    void addBackExplorer(ExplorerI!(T) e){
        activateExplorer(SKeyVal.Master,e.name);
    }
    /// updates a pending operation status, should remove the operation when finished
    void updateEvalStatus(SKey owner,char[] opId,ModifyEvalOp!(T) op,EvalOp!(T).Status newStat){
        if (hasKey(owner)){
            pendingEvals.updateEvalStatus(opId,op,newStat);
        }
    }
    void addEvalOp(SKey s,EvalOp!(T)op){
        if (hasKey(s)){
            localOps.pushBack(op);
        } else {
            silosForKey(s).addEvalOp(s,op);
        }
    }
    /// this should be called by the master process sending it to SKey.All
    /// to receive a new operation to do
    void prepareNextOp(SKey s, int tag){
        if (s==SKeyVal.All){
            foreach (sK,sL;silos){
                sL.prepareNextOp(sK,tag);
            }
        } else if (hasKey(s)){
            ExplorerI!(T) e;
            loop:for (int i=0;i<input.maxExplorationTries;++i){
                char[] explName;
                if (paraEnv.myRank==0 && (!activeExplorers.popFront(explName))){
                    for (int j=0;j<10;++j){
                        waitExplorers.wait();// to do: should have a timeout
                        if (activeExplorers.popFront(explName)) break;
                    }
                }
                mpiBcastT(paraEnv,explName,0,tag);
                if (explName.length>0){
                    synchronized(explorers){
                        e=explorers[explName];
                    }
                    auto res=e.nextOp(&this.addBackExplorer,tag);
                    switch(res){
                    case ExplorerI!(T).ReturnFlag.NoOp:
                        // remove explorer
                        activeExplorers.filterInPlace(/+scope+/ delegate bool(cstring s){
                            return s==explName;
                        });
                        paraEnv.barrier();
                        explorers.removeKey(explName); // this forces explorers to return null only when there are no pending ops...
                        break loop;
                    case ExplorerI!(T).ReturnFlag.SkipOp:
                        if (paraEnv.myRank==0){
                            auto emptyOp=new EvalOp!(T)();
                            registerPendingOp(emptyOp);
                            localOps.pushBack(emptyOp);
                        }
                        break loop;
                    case ExplorerI!(T).ReturnFlag.LocalOp:
                        assert(paraEnv.myRank!=0);
                        if (paraEnv.myRank==0){
                            atomicAdd(pendingEvals.nPending,1);
                        }
                        break loop;
                    default:
                        assert(0);
                    }
                }
                if (explorers.length==0) break loop;
            }
            if (explorers.length==0){
                increaseRunLevel(SKeyVal.Any,RunLevel.WaitPending);
            }
        } else {
            silosForKey(s).prepareNextOp(s,tag);
        }
    }
    /// this should be sent only to the master, and will return a new operation to do
    /// returns an EmptyOp if there is no work
    EvalOp!(T) getNextOp(SKey s){
        if (hasKey(s)){
            assert(hasKey(SKeyVal.Master),"works only on master");
            assert((cast(int)paraEnv.maxTagMask)>0);
            EvalOp!(T) res;
            while (!localOps.popFront(res)){
                int tag;
                do {
                    tag=cast(int)(nextTag.next & paraEnv.maxTagMask);
                } while (tag==0)
                synchronized(this){
                    if (runLevel>=RunLevel.Running) return null;
                    if (nEvals>=input.explorationSteps) {
                        runLevel=RunLevel.WaitPending;
                        waitPendingEnd.checkCondition();
                        return null;
                    }
                    ++nEvals;
                }
            }
            return res;
        } else {
            return silosForKey(s).getNextOp(s);
        }
    }
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(SKey s,Point p){
        assert(0,"to do");
    }
    
    /// returns the key of the silos with the lowest expected load
    SKey nextFreeSilos(SKey s){
        if (hasKey(s)){
            auto w=loads.popWait();
            auto lCtx=this;
            mixin(mkActionMixin("updLoad","w|lCtx|loads",`
            w.load=lCtx.load(w.silosKey);
            loads.push(w);`));
            if (hasKey(w.silosKey)){ // local update
                updLoad();
            } else {
                Task("updLoad",updLoad).autorelease.submitYield();
            }
            return w.silosKey;
        } else {
            return silosForKey(s).nextFreeSilos(s);
        }
    }
    
    int nextSilosRank(int tag){
        ubyte[256] buf;
        auto lMem=LocalMem(buf);
        double myLoad=load(SKeyVal.Any);
        double[] resLoad=lMem.allocArr!(double)(paraEnv.dim);
        mpiAllGatherT(paraEnv,myLoad,resLoad,tag);
        double minLoad=double.min;
        int iVal=0;
        foreach(i,lAtt;resLoad){
            if (minLoad>lAtt){
                minLoad=lAtt;
                iVal=i;
            }
        }
        lMem.deallocArr(resLoad);
        return iVal;
    }
    
    void stop(){
        synchronized(this){
            explorers.clear; // drops all explorers
        }
    }
    
    MainPoint hydrate(V)(DriedPoint!(V) driedPoint){
        auto pSys=refPos.dup();
        assert(driedPoint.pos!is null,"position has to be valid");
        pSys[]=driedPoint.pos;
        DynPVect!(T,DualDxType) mDir;
        if (driedPoint.minDir.isNonNull()){
            mDir=pSys.dynVars.dVarStruct.emptyDualDx();
            mDir[]=minDir;
        }
        MainPoint res=new MainPoint(this,driedPoint.point,pSys,driedPoint.gFlags|GFlags.LocalCopy,
            driedPoint.explorationSize,mDir,driedPoint.minDirScale,driedPoint.exploredDirs,driedPoint.neighbors);
        return res;
    }
    
    // ExplorationObserverI
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy){
        if (hasKey(s)) {
            auto mp=localPoints[p];
            if (! mp.isLocalCopy){
                throw new Exception("Local copy in PNetSilos.addEnergyEvalLocal",__FILE__,__LINE__);
            }
            mp.addEnergyEvalLocal(energy);
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.addEnergyEvalLocal(s,p,energy);
            });
        } else {
            silosForKey(s).addEnergyEvalLocal(s,p,energy);
        }
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN does not replace the previous value
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){
        if (s==SKeyVal.All){
            foreach (sK,sL;silos){
                sL.addGradEvalLocal(sK,p,pSys);
            }
        } else if (hasKey(s)) {
            auto mp=localPoints[p];
            if (mp.isLocalCopy){
                throw new Exception("Local copy in PNetSilos.addGradEvalLocal",__FILE__,__LINE__);
            }
            synchronized(mp){
                auto e=mp.pos.dynVars.potentialEnergy;
                mp.pos[]=pSys;
                if (isNaN(mp.pos.dynVars.potentialEnergy)) mp.pos.dynVars.potentialEnergy=e;
            }
            mp.didGradEval();
        } else {
            silosForKey(s).addGradEvalLocal(s,p,pSys);
        }
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void publishPoint(SKey s,SKey owner,Point newPoint,PSysWriter!(T) pos,T pSize,uint flags){
        if (s==SKeyVal.All){
            foreach (sK,sL;silos){
                sL.publishPoint(sK,owner,newPoint,pos,pSize,flags);
            }
        } else if (hasKey(s)) {
            bool collision=false;
            synchronized(this.owner){
                nextPntNr.ensure(newPoint.data>>12);
                auto oldV=newPoint in this.owner;
                if (oldV !is null && (*oldV)!=owner){
                    collision=true;
                } else {
                    this.owner[newPoint]=owner;
                }
            }
            if (collision){
                // doubly owned point (collision): drop
                publishCollision(SKeyVal.All,newPoint);
                return;
            }
            silosForKey(owner).didLocalPublish(owner,newPoint,s);
            foreach(p,mp;localPoints){
                auto newPos=refPos.dynVars.dVarStruct.emptyX();
                newPos[]=pos.x;
                if (mp!is null){
                    mp.checkIfNeighbor(newPos,newPoint,pSize);
                }
                newPos.giveBack();
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.publishPoint(s,owner,newPoint,pos,pSize,flags);
            });
        } else {
            silosForKey(s).publishPoint(s,owner,newPoint,pos,pSize,flags);
        }
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to silos
    void neighborHasEnergy(SKey s,Point p,Point[] neighbors,Real energy){
        if (hasKey(s)){
            foreach (nP;neighbors){
                auto mp=localPoints[nP];
                mp.addEnergyEvalOther(p,energy);
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.neighborHasEnergy(s,p,neighbors,energy);
            });
        } else {
            silosForKey(s).neighborHasEnergy(s,p,neighbors,energy);
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors might be restricted to silos or not
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, Real energy){
        if (hasKey(s)){
            foreach (nP;neighbors){
                auto mp=localPoints[nP];
                mp.addGradEvalOther(p,energy);
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.neighborHasGradient(s,p,neighbors,energy);
            });
        } else {
            silosForKey(s).neighborHasGradient(s,p,neighbors,energy);
        }
    }
    
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPoint(SKey s,Point point,SKey owner){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.finishedExploringPoint(sK,point,owner);
            }
        } else if (hasKey(s)){
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.finishedExploringPoint(s,point,owner);
            });
        } else {
            silosForKey(s).finishedExploringPoint(s,point,owner);
        }
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey s,Point point){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.publishCollision(sK,point);
            }
        } else if (hasKey(s)){
            MainPointI!(T) mp;
            synchronized(localPoints){
                auto mpp=point in localPoints;
                if (mpp !is null)
                    mp=*mpp;
            }
            if (mp!is null){
                PointBcastProgress!(T) pp;
                synchronized(publishingPoints){
                    pp=publishingPoints[point];
                }
                pp.communicateCollision();
                mp.drop(); // redundant...
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.publishCollision(s,point);
            });
        } else {
            silosForKey(s).publishCollision(s,point);
        }
    }
    /// informs that source has processed point p0
    void didLocalPublish(SKey s,Point p0,SKey source){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.didLocalPublish(sK,p0,source);
            }
        } else if (hasKey(s)){
            PointBcastProgress!(T) pp;
            synchronized(publishingPoints){
                pp=publishingPoints[p0];
            }
            pp.communicatePublish(source);
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.didLocalPublish(s,p0,source);
            });
        } else {
            silosForKey(s).didLocalPublish(s,p0,source);
        }
    }
    /// increases the runLevel on all silos, i.e. you should call it only with SKeyVal.All
    /// (at the moment there is no support for dynamic adding/removal of silos)
    void increaseRunLevel(SKey s,RunLevel newRunLevel){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.increaseRunLevel(sK,newRunLevel);
            }
        } else if (hasKey(s)){
            bool runLevelChanged=false;
            synchronized(this){
                if (runLevel<newRunLevel){
                    runLevel=newRunLevel;
                    runLevelChanged=true;
                }
            }
            if (runLevelChanged){ // square b-cast, change??
                waitPendingEnd.checkCondition();
                notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                    obs.increaseRunLevel(s,newRunLevel);
                });
            }
        } else {
            assert(0,"should be called only with SKeyVal.All");
        }
    }

    // PNetSilosI

    /// load (usage) of the silos in some arbitrary units
    Real load(SKey s){
        if (hasKey(s)){
            return 0.001*localPoints.length;
        } else {
            return silosForKey(s).load(s);
        }
    }
    
    /// energy for the local points (NAN if not yet known)
    Real[] energyForPointsLocal(SKey s,Point[] points,Real[] ens){
        if (hasKey(s)){
            if (ens.length<points.length)
                ens.length=points.length;
            foreach(i,p;points){
                auto mp=localPoints[p];
                ens[i]=mp.energy;
            }
            return ens;
        } else {
            return silosForKey(s).energyForPointsLocal(s,points,ens);
        }
    }
    /// energy for the points (NAN if not yet known)
    Real[] energyForPoints(SKey s,Point[]points,Real[] ens){
        if (hasKey(s)){
            if (ens.length<points.length)
                ens.length=points.length;
            auto pO=new PointToOwner(&this.ownerOfPoint,delegate Point(size_t i){
                return points[i];
            }, points.length,false);
            Real[128] buf;
            auto lMem=LocalMem(buf);
            foreach (iter;pO){
                auto enLoc=lMem.allocArr!(Real)(iter.localUb-iter.localLb);
                enLoc=energyForPointsLocal(iter.owner,iter.points,enLoc);
                foreach(i,j;iter.idx){
                    ens[j]=enLoc[i];
                }
                lMem.deallocArr(enLoc);
            }
            return ens;
        } else {
            return silosForKey(s).energyForPoints(s,points,ens);
        }
    }
    /// returns a snapshot of the given point
    DriedPoint!(T)mainPoint(SKey s,Point p){
        if (hasKey(s)){
            auto ownr=ownerOfPoint(p);
            if (hasKey(ownr)){
                return mainPointLocal(s,p);
            }
            return silosForKey(ownr).mainPointLocal(ownr,p);
        } else {
            return silosForKey(s).mainPoint(s,p); // this allows access to non public points...
        }
    }
    /// returns a snapshot of the given point that is local to this silos
    DriedPoint!(T)mainPointLocal(SKey s,Point p){
        if (hasKey(s)){
            auto mp=localPoints[p];
            return mp.driedPoint;
        } else {
            silosForKey(s).mainPointLocal(s,p);
        }
    }
    /// tells the local points neighs that they have p0 as neighbor
    void addPointToLocalNeighs(SKey s,Point p0,Point[]neighs){
        if (hasKey(s)){
            foreach(np;neighs){
                auto mp=localPoints[np];
                mp.addNeighbor(p0);
            }
        } else {
            silosForKey(s).addPointToLocalNeighs(s,p0,neighs);
        }
    }
    /// tells the local point p0 that neighDirs (computed when p0 gradient was hadGrad) should be added to it
    bool addNeighDirsToLocalPoint(SKey s,Point p0,PointAndDir[]neighDirs,DirDistances!(T)[]dirDists,bool hadGrad){
        if (hasKey(s)){
            auto mp=localPoints[p0];
            return mp.addNeighbors(neighDirs,dirDists,hadGrad);
        } else {
            silosForKey(s).addNeighDirsToLocalPoint(s,p0,neighDirs,dirDists,hadGrad);
        }
    }
    /// creates a new point located at newPos in this silos, the point is not yet broadcasted
    /// not all silos might support creation of local points, use nextFreeSilos to get a silos
    /// thas supports it, use a RemoteSilosOpI and executeLocally to create, setup & bcast it 
    MainPointI!(T) newPointAt(DynPVector!(T,XType) newPos){
        MainPoint!(T) newP=pPool.getObj();
        newP.pos.checkX();
        newP.pos.dynVars.x[]=newPos;
        synchronized(localPoints){ // don't allow null mainpoints...
            Point np;
            synchronized(owner){ // generate a point with a small propability of collisions
                while(1){
                    ushort r;
                    rand()(r);
                    np=Point(((nextPntNr.next())<<12)|cast(ulong)(r&0xFFF));
                    localPointsKeys.appendEl(np);
                    auto oldV=np in owner;
                    if (oldV is null){
                        owner[np]=this._key;
                        break;
                    }// else retry immediately in case of collision
                }
            }
            newP._point=np;
            localPoints[np]=newP;
        }
        return newP;
    }
    /// operation to be executed on the given silos
    void executeLocally(SKey s,RemoteSilosOpI!(T) op){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.executeLocally(sK,op);
            }
        } else if (hasKey(s)){
            op.workOn(this);
        } else {
            silosForKey(s).executeLocally(s,op);
        }
    }
    
    // LocalSilosI
    /// adds an extra observer that will be notified about the network changes
    void addObserver(ExplorationObserverI!(T) o){
        observers.append(o);
    }
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o){
        observers.filterInPlace(delegate bool(ExplorationObserverI!(T) i){ return i!is o; });
    }
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o){
        synchronized(explorers){
            explorers[o.name]=o;
        }
    }
    /// removes the given explorer
    void rmExplorerNamed(char[] name){
        synchronized(explorers){
            explorers.removeKey(name);
        }
    }
    void notifyLocalObservers(void delegate(ExplorationObserverI!(T))op){
        foreach(obs;observers){
            op(obs);
        }
    }
    /// linear communicator (valid only inside the real silos, not in the clients)
    LinearComm paraEnv(){
        return silosParaEnv;
    }
    SKey pointOwner(SKey s,Point p){
        if (hasKey(s)){
            synchronized(owner){
                return owner[p];
            }
        } else {
            silosForKey(s).pointOwner(s,p); // can be useful for not yet broadcasted points
        }
    }
    SKey ownerOfPoint(Point p){
        return pointOwner(SKeyVal.Any,p);
    }
    /// internal method returning the silos for the given key
    protected PNetSilosI!(T) silosForKey(SKey k){
        if (k==SKeyVal.Invalid||k==SKeyVal.All){
            throw new Exception(((k==SKeyVal.Invalid)?"invalid SKey"[]:"broadcast key not acceptable in this context"[]),__FILE__,__LINE__);
        }
        if (k==SKeyVal.Master) k=sKeys[0];
        if (hasKey(k)) return this;
        return silos[k];
    }
    /// local point mainpoint (the real reference point)
    MainPointI!(T) mainPointL(Point p){
        return localPoints[p];
    }
    /// a local point (possibly a copy), is retained, and needs to be released (thus the create in the name)
    /// the time t is used to decide if a cached version can be used
    MainPointI!(T)createLocalPoint(Point p,ev_tstamp t){
        if (hasKey(ownerOfPoint(p))){
            auto mp=localPoints[p];
            mp.retain;
            return mp;
        }
        CachedPoint!(T) cachedP;
        CachedPoint!(T) * pC;
        synchronized(localCache){
            pC= p in localCache;
            if (pC!is null){
                cachedP=*pC;
            }
        }
        if (pC is null){
            cachedP.mainPoint=pPool.getObj();
            cachedP.mainPoint.pos.checkX();
            cachedP.mainPoint._point=p;
            cachedP.mainPoint._gFlags|=GFlags.LocalCopy;
        } else if (cachedP.lastSync>t){
            cachedP.mainPoint.retain;
            return cachedP.mainPoint;
        }
        cachedP.lastSync=ev_time();
        cachedP.mainPoint[]=mainPoint(ownerOfPoint(p),p);
        synchronized(localCache){
            localCache[p]=cachedP;
        }
        return cachedP.mainPoint;
    }
    /// makes a point "public" informing other silos that that region has been explored
    /// returns the published point (that might be different if there was a collision)
    MainPointI!(T) bcastPoint(MainPointI!(T) p){
        auto bcastP=PointBcastProgress!(T)(p);
        bcastP.publish();
        return bcastP.finalName();
    }
    /// drops a cached point (the point is not in use anymore)
    void dropCachedPoint(MainPointI!(T)p){
        synchronized(localCache){
            localCache.removeKey(p.point);
        }
    }
    /// registers a pending operation (sets id,...)
    void registerPendingOp(EvalOp!(T)op){
        pendingEvals.addPendingOp(op);
    }
    mixin(rpcMixin("dchem.PNetSilos!("~T.stringof~")", "PNetSilosI!("~T.stringof~")",silosMethodsStr));
    DefaultVendor mainVendor;
    char[] silosCoreUrl(){
        return mainVendor.proxyObjUrl();
    }
    this(SilosGen sGen,char[] nameId,SKey key,ParticleSys!(T) refPos=null,ConstraintI!(T) constraints=null,
        CharSink log=sout.call,NotificationCenter _nCenter=null,PoolI!(MainPoint!(T)) pPool=null){
        this.input=sGen;
        assert(this.input!is null);
        this.runLevel=RunLevel.Setup;
        this.nextLocalId=UniqueNumber!(ulong)(1);
        this.nextPntNr=UniqueNumber!(ulong)(3);
        this.nameId=nameId;
        this._key=key;
        this.silos=new HashMap!(SKey,PNetSilosI!(T))();
        this.silos[this._key]=this;
        this.silosParaEnv=mpiWorld;
        this.loads=new MinHeapSync!(LoadStats)();
        this.localPointsKeys=new BatchedGrowableArray!(Point,batchSize)();
        this.owner=new HashMap!(Point,SKey)();
        this.localCache=new HashMap!(Point,CachedPoint!(T))();
        this.pendingEvals=new PendingEvals!(T)(this);
        this.localOps=new Deque!(EvalOp!(T))();
        this.nPendingOps=0;
        this.localPoints=new DenseLocalPointArray!(MainPointI!(T))(); /// local points (position,...)
        this._rand=new RandomSync();
        this._refPos=refPos;
        this._constraints=constraints;
        this.log=log;
        this._nCenter=nCenter;
        if (nCenter is null) this._nCenter=new NotificationCenter();
        this.explorers=new HashMap!(char[],ExplorerI!(T))();
        this.observers=new Deque!(ExplorationObserverI!(T))();
        this.evaluator=cast(Method)this.input.evaluator.contentObj;
        if (this.input.evaluatorTask!is null) {
            this.evaluatorTask=cast(RemoteCCTask)this.input.evaluatorTask.contentObj;
        }
        if (this.evaluatorTask is null){
            // add this dependency??
            // this.evaluatorTask=new WorkAsker();
            // this.evaluatorTask.precision=this.input.precision;
        }
        this.pPool=pPool;
        if (pPool is null) this.pPool=cachedPool(&this.allocPoint);
        foreach(e;input.explorers){
            auto o=e.contentObj;
            if (cast(ExplorerGen)o !is null){
                auto expl=explorerT!(T)(cast(ExplorerGen)o,this);
                explorers[expl.name]=expl;
                observers.append(expl);
            } else {
                observers.append(observerT!(T)(cast(ExplorationObserverGen)o,this));
            }
        }
        this.activeExplorers=new Deque!(string)();
        foreach (k,e;explorers){
            this.activeExplorers.pushBack(k);
        }
        this.waitExplorers=new WaitCondition(&this.hasActiveExpl);
        this.waitPendingEnd=new WaitCondition(&this.isAtPendingEnd);
        this.loaders=new Deque!(SilosWorkerI!(T))();
        this.monitors=new Deque!(SilosWorkerI!(T))();
        this.finishers=new Deque!(SilosWorkerI!(T))();
        this.mainVendor=new DefaultVendor(this);
        ProtocolHandler.defaultProtocol.publisher.publishObject(mainVendor,"silos",true);
        foreach(l;input.loaders){
            auto swGen=cast(SilosWorkerGen)(l.contentObj());
            auto sw=silosWorkerT!(T)(swGen);
            loaders.append(sw);
        }
        foreach(m;input.monitors){
            auto swGen=cast(SilosWorkerGen)(m.contentObj());
            auto sw=silosWorkerT!(T)(swGen);
            monitors.append(sw);
        }
        foreach(f;input.finishers){
            auto swGen=cast(SilosWorkerGen)(f.contentObj());
            auto sw=silosWorkerT!(T)(swGen);
            finishers.append(sw);
        }
    }
}

PNetSilos!(Real) dummyPNetSilos;
