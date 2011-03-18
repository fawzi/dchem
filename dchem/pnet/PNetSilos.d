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
import dchem.pnet.WaitOp;
import blip.container.AtomicSLink;
import blip.util.RefCount;
import blip.core.Array:sort;
//import dchem.pnet.WorkAsker;

/// help structure 
struct CachedPoint(T){
    MainPoint!(T) mainPoint;
    ev_tstamp lastSync;
}

struct LoadStats{
    Real load;
    SKey silosKey;
    mixin(serializeSome("LoadStats","statistic about the load (usage) of a silos","load|silosKey"));
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

class SilosGen:Sampler {
    PNetSilos!(Real)[] silosReal;
    PNetSilos!(LowP)[] silosLowP;
    
    /// number of exploration steps to perform
    long explorationSteps=long.max;
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
    /// eagerness in calculating the gradients
    char[] gradEagerness="OnRequest";
    /// list of observers and explorers
    InputField[] explorers;
    /// list of loaders (called after the explorer are set up, before exploring)
    InputField[] loaders; // use simply SilosWorkerGen[]??
    /// list of monitors, called at regular intervals during the exploration
    InputField[] monitors; // use simply SilosWorkerGen[]??
    /// list of finishers, called after the exploration is finished
    InputField[] finishers; // use simply SilosWorkerGen[]??
    
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
    
    /// the time to wait when no task is available
    double noTaskWaitTime=10;
    
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
    mixin(serializeSome("dchem.Silos",`A silos, the place where all points are stored and observers and explorers live.`,
        `noTaskWaitTime: the number of seconds to wait when no task is available (10)
        evaluator: the method to perform energy evaluations
        evaluatorTask: the task to execute (normally to start clients that will perform evaluations, if a SilosConnectorI then it is passed silos and accuracy before starting)
        precision: the precision with which the points are stored
        gradEagerness: the eagerness in calculating gradients (OnRequest,Speculative,Always)
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
    mixin printOut!();
    this(){}

    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        if (evaluator is null || cast(Method)evaluator.contentObj is null){
            log("evaluator should be valid and be a Method in field ")(myFieldName)("\n");
            res=false;
        }
        if (evaluatorTask !is null && cast(Sampler)evaluatorTask.contentObj is null){
            log("evaluatorTask if given should be a Sampler in field ")(myFieldName)("\n");
            res=false;
        }
        if (precision!="LowP" && precision!="Real"){
            log("precision should be either LowP or Real, not '")(precision)("' in field ")(myFieldName)("\n");
            res=false;
        }
        if (gradEagerness!="OnRequest" && gradEagerness!="Speculative" && gradEagerness!="Always"){
            log("gradEagerness should be either OnRequest or Speculative or Always, not '")(gradEagerness)("' in field ")(myFieldName)("\n");
            res=false;
        }
        foreach(e;explorers){
            auto o=cast(Object)e.content;
            if (cast(ExplorationObserverGen)o is null){
                log("explorers should be of either explorers or observer, not '")(o.classinfo.name)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(l;loaders){
            auto o=cast(Object)l.content;
            if (cast(SilosWorkerGen)o is null){
                log("loaders should be a SilosWorker not '")(o.classinfo.name)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(m;monitors){
            auto o=cast(Object)m.content;
            if (cast(SilosWorkerGen)o is null){
                log("monitors should be a SilosWorker not '")(o.classinfo.name)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(f;finishers){
            auto o=cast(Object)f.content;
            if (cast(SilosWorkerGen)o is null){
                log("finishers should be a SilosWorker not '")(o.classinfo.name)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        return res;
    }
    
    GradEagerness gEagerness(){
        switch(gradEagerness){
            case "OnRequest":
                return GradEagerness.OnRequest;
            case "Speculative":
                return GradEagerness.Speculative;
            case "Always":
                return GradEagerness.Always;
            default:
                assert(0,"invalid gradEagerness");
        }
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
    long waitingName=1,waitingNeigh=1;
    WaitConditionT!(false) nameKnown;
    PointBcastProgress nextPoint;
    PoolI!(PointBcastProgress) pool;
    this(PoolI!(PointBcastProgress)p){
        pool=p;
        this.nameKnown=new WaitConditionT!(false)(&this.nameIsKnown);
    }
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
        refCount=1;
        waitingName=1;
        waitingNeigh=1;
        this.localContext=cast(PNetSilos!(T))cast(Object)p.localContext;
        assert(this.localContext!is null);
    }
    static PointBcastProgress opCall(MainPointI!(T)p){
        auto nP=gPool.getObj();
        nP.reset(p);
        return nP;
    }
    void publish(){
        retain();
        synchronized(localContext.publishingPoints){
            this.localContext.publishingPoints[pointToBCast.point]=this;
        }
        foreach(sK,sL;localContext.silos){
            synchronized(this){
                ++waitingName;
                ++waitingNeigh;
            }
            sL.publishPoint(sK,localContext._key,pointToBCast.point,pSysWriter(pointToBCast.pos),
                pointToBCast.explorationSize,pointToBCast.gFlags);
        }
        bool last;
        synchronized(this){
            --waitingName;
            --waitingNeigh;
            last= (waitingNeigh==0);
        }
        nameKnown.checkCondition();
        if (last){
            pointToBCast.bcastLevel(1);
            release();
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
            auto np=localContext.newPointAt(pointToBCast.pos.dynVars.x,Point(0));
            np[]=dP;
            nextPP=PointBcastProgress(np);
            nextPP.publish();
        }
        bool last=false;
        synchronized(this){
            if (nextPP!is null) nextPoint=nextPP;
            --waitingName;
            --waitingNeigh;
            if (waitingNeigh==0) last=true;
        }
        if (last){
            pointToBCast.bcastLevel(1);
            release();
        }
        nameKnown.checkCondition();
    }
    void communicatePublish(SKey origin,int level){
        bool last=false;
        switch(level){
        case 0:
            synchronized(this){
                --waitingName;
            }
            assert(nameKnown!is null);
            nameKnown.checkCondition();
            break;
        case 1:
            synchronized(this){
                --waitingNeigh;
                if (waitingNeigh==0) last=true;
            }
            if (waitingNeigh<0){
                throw new Exception("waitingNeigh<0",__FILE__,__LINE__);
            }
            if (last) {
                pointToBCast.bcastLevel(1);
                release();
            }
            break;
        default:
            assert(0);
        }
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
        res.bcastLevel(0);
        localContext.publishedLocalPoint(SKeyVal.Any,res.point);
        release();
        return res;
    }
    void release0(){
        synchronized(localContext.publishingPoints){
            localContext.publishingPoints.removeKey(pointToBCast.point);
        }
        if (pool!is null) pool.giveBack(this);
    }
    mixin RefCountMixin!();
}

class PendingEvals(T) {
    PNetSilos!(T) silos;
    HashMap!(char[],EvalOp!(T)) pendingOps;
    UniqueNumber!(ulong) nextLocalId;
    
    this(PNetSilos!(T) silos){
        this.silos=silos;
        this.pendingOps=new HashMap!(char[],EvalOp!(T))();
        this.nextLocalId=UniqueNumber!(ulong)(1);
    }
    size_t nPending(){
        synchronized(this){
            return pendingOps.length;
        }
    }
    void addPendingOp(EvalOp!(T) op){
        auto id=collectAppender(delegate void(CharSink s){
            dumper(s)(silos._key)("-")(nextLocalId.next());
        });
        version (TrackPNet) silos.logMsg(delegate void(CharSink s){
            dumper(s)("addPendingOp(")(op)(") with id:")(id);
        });
        op.initMaster(id,silos,silos._key);
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
                silos.addEvalOp(SKeyVal.Master,opAtt,false);
                resubmitted=true;
            }
        }
        version (TrackPNet) silos.logMsg(delegate void(CharSink s){
            dumper(s)("updateStatus(")(opId)(",")((op is null)?"*null*":"op")(",")(newStat)(")");
        });
        switch (newStat){
            case EvalOp!(T).Status.ToDo:
                assert(0);
            case EvalOp!(T).Status.InProgress:
                break;
            case EvalOp!(T).Status.Failure:
            case EvalOp!(T).Status.Success:
                if (!resubmitted){
                    if (opAtt){
                        synchronized(pendingOps){
                            pendingOps.removeKey(opId);
                        }
                        if (silos.paraEnv.myRank!=0){
                            silos.updateEvalStatus(SKeyVal.Master,opId,null,newStat);
                        }
                    }
                    if (silos.paraEnv.myRank==0){
                        atomicAdd(silos.nPendingOps,-1);
                        silos.waitPendingEnd.checkCondition();
                    }
                }
                break;
            default:
                assert(0);
        }
    }
}

final class PNetSilos(T): LocalSilosI!(T){
    RunLevel runLevel;
    SilosGen input;
    UniqueNumber!(ulong) nextLocalId;
    UniqueNumber!(ulong) nextPntNr;
    UniqueNumber!(uint)  nextTag;
    long nPendingOps=0; /// global number of pending ops (on the master only)
    string nameId;
    SKey _key;
    bool hasKey(SKey k){ return _key==k||k==SKeyVal.Any||(k==SKeyVal.Master && paraEnv.myRank==0); }
    HashMap!(SKey,PNetSilosI!(T)) silos; // contains self
    MinHeapSync!(LoadStats) loads; 
    BatchedGrowableArray!(Point,batchSize) _localPointsKeys; // indexes for local owned points
    BatchedGrowableArray!(Point,batchSize) localPointsKeys(){
        return _localPointsKeys;
    }
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
    HashMap!(char[],ExplorerI!(T)) explorers; /// all explorers (so that concurrent nextOp execution doen't give problems). If you have many you should probably improve this
    HashMap!(char[],ExplorerI!(T)) generatingExplorers; /// the explorers that might generate new points (but might be suspended)
    Deque!(char[]) activeExplorers; /// explorers that can generate new points
    WaitCondition waitExplorers;
    WaitCondition waitPendingEnd;
    Deque!(EvalOp!(T)) localOps;
    /// parallel enviroment of silos
    LinearComm silosParaEnv;
    /// observers
    alias SLinkT!(ExplorationObserverI!(T)) ObserversEl;
    ObserversEl* observers; // could use a better persistent structure... but it isn't critical
    HashMap!(string,SilosComponentI!(T)) _components; // ideally should also be persistent...
    
    /// evaluator
    Method evaluator;
    Sampler evaluatorTask;
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
        return activeExplorers.length>0 || generatingExplorers.length==0;
    }
    bool isStopping(){
        return runLevel>RunLevel.WaitPending;
    }
    /// writes out a log message
    void logMsg(void delegate(CharSink)writer){
        char[512] buf;
        auto gArr=lGrowableArray(buf,0);
        auto s=dumper(&gArr.appendArr);
        s("<PNetSilosLog id=\"")(_key)("\" time=\"")(ev_time())("\" task=\"")(taskAtt.val)("\">\n  ");
        indentWriter("  ",s.call,writer);
        s("\n</PNetSilosLog>\n");
        log(gArr.data);
        gArr.deallocData;
    }
    /// writes out a log message
    void logMsg1(char[]msg){
        logMsg(delegate void(CharSink s){ s(msg); });
    }
    /// constraints of the current system
    ConstraintI!(T) constraints(){
        return _constraints;
    }
    /// how eagerly the gradient is calculated
    GradEagerness gradEagerness(){
        return input.gEagerness(); // cache?
    }
    /// if the gradient should be speculatively calculated. calls addEnergyEvalLocal (in background)
    /// (should be called on the owner of the point only)
    bool speculativeGradientLocal(SKey s,Point p,Real energy,Real energyError){
        if (hasKey(s)){
            bool res=false;
            int i=0;
            auto pAtt=mainPointL(p);
            if (pAtt.setEnergy(energy,energyError)){
                Task("didEnergyEval",&pAtt.didEnergyEval).autorelease.submit(defaultTask);
            }
            if ((pAtt.gFlags&GFlags.GradientInfo)==0){ // not already in progress or done...
                notifyLocalObservers(delegate void(ExplorationObserverI!(T) obs){
                    if (obs.speculativeGradientLocal(_key,p,energy,energyError)) res=true;
                });
            }
            if (res) res=pAtt.shouldCalculateGradient();
            return res;
        } else {
            return silosForKey(s).speculativeGradientLocal(s,p,energy,energyError);
        }
    }
    
    mixin(realFromInput(propertiesList));
    
    Real[char[]] propertiesDict(SKey s){
        if (hasKey(s)){
            auto res=propertiesDict0();
            res["gradEagerness"]=cast(Real)cast(uint)gradEagerness();
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
            volatile bool hasOps= (paraEnv.myRank==0 && nPendingOps==0);
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
        assert((!res)||pendingEvals.nPending==0);
        return res;
    }
    mixin(descSome("dchem.PNetSilos_"~T.stringof,`_key|runLevel|nextLocalId|nextPntNr|nameId`));
    mixin printOut!();
    
    void startAskers(){
        if (evaluatorTask is null) return;
        auto silosC=cast(SilosConnectorI)cast(Object)evaluatorTask;
        if (silosC!is null){
            silosC.setConnectionAndPrecision(silosCoreUrl,input.precision);
        }
        logMsg1("starting evaluatorTask");
        evaluatorTask.run(silosParaEnv,log);
        logMsg1("evaluatorTask finished");
    }
    
    void start(LinearComm pEnv, CharSink log){
        silosParaEnv=pEnv;
        pEnv.barrier();
        logMsg1("silos check SKeys");
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
        pEnv.barrier();
        logMsg1("silos setting up evaluator");
        evaluator.setup(pEnv,log);
        CalculationContext cInstance;
        if (_refPos is null || _constraints is null){
            cInstance=evaluator.getCalculator(true,[]);
        }
        logMsg1("silos setting up refPos");
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
        logMsg1("silos setting up constraints");
        if (_constraints is null){
            auto cGen=cInstance.constraintGen;
            if (cGen is null) {
                _constraints=new NoConstraint!(T)();
            } else {
                _constraints=constraintT!(T)(cGen,_refPos);
            }
        }
        _constraints.applyR(_refPos);
        if (cInstance!is null){
            cInstance.giveBack();
        }
        logMsg1("silos will run the loaders");
        // run the loaders
        foreach(sw;loaders){
            sw.workOn(this);
        }
        paraEnv.barrier();
        logMsg1("silos did run the loaders");
        runLevel=RunLevel.Running;
        int nOps=localOps.length;
        paraEnv.bcast(nOps,0);
        if (generatingExplorers.length==0 && nOps==0){
            if (paraEnv.myRank==0) log("no explorers, and no pendingOps, stopping (to avoid this add a dchem.WaitExplorer)\n");
            increaseRunLevel(SKeyVal.All,RunLevel.WaitPending);
        }
        /// start the askers
        Task("startAskers",&startAskers).autorelease.submitYield();
    }
    void run(LinearComm pEnv, CharSink log){
        logMsg(delegate void(CharSink s){
            dumper(s)("PNetSilos!(")(T.stringof)(") was generated from field '")(input.myFieldName)("'");
        });
        start(pEnv,log);
        logMsg1("silos started waiting");
        waitPendingEnd.wait();
        logMsg1("silos finished waiting");
        if (paraEnv.myRank==0){
            increaseRunLevel(SKeyVal.All,RunLevel.Stopping);
        }
        pEnv.barrier();
        synchronized(pendingEvals){
            if (runLevel>RunLevel.Stopping) return;
            runLevel=RunLevel.Stopping;
        }
        logMsg1("silos will run the finishers");
        // run the finishers
        foreach(sw;finishers){
            sw.workOn(this);
        }
        // stop components
        foreach (comp;components.dup){
            comp.stop();
        }
        pEnv.barrier();
        logMsg1("silos did run the finishers");
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
    void addEvalOp(SKey s,EvalOp!(T)op,bool incrementNPending){
        if (hasKey(s)){
            if (op.id.length==0){
                throw new Exception("non registred operation in addEvalOp",__FILE__,__LINE__);
            }
            localOps.pushBack(op);
            if (incrementNPending) atomicAdd(nPendingOps,1);
        } else {
            silosForKey(s).addEvalOp(s,op,incrementNPending);
        }
    }
    /// this should be called by the master process sending it to SKey.All
    /// to receive a new operation to do
    void prepareNextOp(SKey s, int tag){
        if (s==SKeyVal.All){
            SKey myK;
            foreach (sK,sL;silos){
                if (!hasKey(sK)){
                    sL.prepareNextOp(sK,tag);
                } else {
                    myK=sK;
                }
            }
            assert(hasKey(myK));
            prepareNextOp(myK,tag);
        } else if (hasKey(s)){
            ExplorerI!(T) e;
            loop:for (int i=0;i<input.maxExplorationTries;++i){
                char[] explName;
                if (paraEnv.myRank==0){
                    if (!activeExplorers.popFront(explName)){
                        for (int j=0;j<10;++j){
                            waitExplorers.wait();// to do: should have a timeout
                            if (activeExplorers.popFront(explName)) break;
                        }
                    }
                    explName=explName.dup; // bcast needs a mutable copy...
                }
                mpiBcastT(paraEnv,explName,0,tag);
                logMsg(delegate void(CharSink s){
                    dumper(s)("prepareNextOp choosenExplorer:")(explName);
                });
                if (explName.length>0){
                    synchronized(explorers){
                        e=explorers[explName];
                    }
                    auto res=e.nextOp(&this.addBackExplorer,tag);
                    switch(res){
                    case ExplorerI!(T).ReturnFlag.NoOp:
                        version(TrackPNet) logMsg(delegate void(CharSink s){
                            dumper(s)("removing explorer");
                        });
                        // remove explorer
                        activeExplorers.filterInPlace(/+scope+/ delegate bool(cstring s){
                            return s==explName;
                        });
                        paraEnv.barrier();
                        synchronized(generatingExplorers){
                            generatingExplorers.removeKey(explName); // this forces explorers to return null only when there are no pending ops...
                        }
                        break loop;
                    case ExplorerI!(T).ReturnFlag.SkipOp:
                        version(TrackPNet) logMsg(delegate void(CharSink s){
                            dumper(s)("skipping explorer");
                        });
                        if (paraEnv.myRank==0){
                            auto emptyOp=new WaitOp!(T)(input.noTaskWaitTime);
                            registerPendingOp(emptyOp);
                            addEvalOp(SKeyVal.Master,emptyOp,true);
                        }
                        break loop;
                    case ExplorerI!(T).ReturnFlag.LocalOp:
                        version(TrackPNet) logMsg(delegate void(CharSink s){
                            dumper(s)("added op");
                        });
                        if (paraEnv.myRank==0){
                            atomicAdd(nPendingOps,1);
                        }
                        break loop;
                    default:
                        assert(0);
                    }
                }
                if (generatingExplorers.length==0) break loop;
            }
            if (generatingExplorers.length==0){
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
            version(TrackPNet) logMsg1("starting getNextOp");
            assert(hasKey(SKeyVal.Master),"works only on master");
            assert((cast(int)paraEnv.maxTagMask)>0);
            EvalOp!(T) res;
            while (!localOps.popFront(res)){
                int tag;
                do {
                    tag=cast(int)(nextTag.next & paraEnv.maxTagMask);
                } while (tag==0)
                synchronized(this){
                    if (runLevel>RunLevel.Running) {
                        version(TrackPNet) logMsg(delegate void(CharSink s){ dumper(s)("getNextOp runLevel>running, will return null"); });
                        return null;
                    }
                    if (nEvals>=input.explorationSteps) {
                        runLevel=RunLevel.WaitPending;
                        waitPendingEnd.checkCondition();
                        version(TrackPNet) logMsg(delegate void(CharSink s){ dumper(s)("getNextOp will return null"); });
                        return null;
                    }
                    ++nEvals;
                }
                prepareNextOp(SKeyVal.All,tag);
            }
            version(TrackPNet) logMsg(delegate void(CharSink s){ dumper(s)("getNextOp will return ")(res); });
            return res;
        } else {
            return silosForKey(s).getNextOp(s);
        }
    }
    /// called when an evaluation fails
    void evaluationFailed(SKey s,Point p){
        assert(0,"to do");
    }
    /// checks if the local point is somehow invalid and should better be skipped
    bool shouldFilterLocalPoint(SKey s,Point p){
        if (hasKey(s)){
            bool res=false;
            notifyLocalObservers(delegate void(ExplorationObserverI!(T) obs){
                if (obs.shouldFilterLocalPoint(_key,p)) res=true;
            });
            return res;
        } else {
            return silosForKey(s).shouldFilterLocalPoint(s,p);
        }
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
            generatingExplorers.clear; // drops all explorers
        }
    }
    
    MainPoint hydrate(V)(DriedPoint!(V) driedPoint){
        auto pSys=refPos.dup();
        assert(driedPoint.pos!is null,"position has to be valid");
        pSys[]=driedPoint.pos;
        DynPVect!(T,DualDxType) mDir;
        if (!driedPoint.minDir.isDummy()){
            mDir=pSys.dynVars.dVarStruct.emptyDualDx();
            mDir[]=minDir;
        }
        MainPoint res=new MainPoint(this,driedPoint.point,pSys,driedPoint.gFlags|GFlags.LocalCopy,
            driedPoint.explorationSize,mDir,driedPoint.minDirScale,driedPoint.exploredDirs,driedPoint.neighbors);
        return res;
    }
    
    // ExplorationObserverI
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy,Real energyError){
        if (hasKey(s)) {
            auto mp=localPoints[p];
            if (mp.isLocalCopy){
                throw new Exception("Local copy in PNetSilos.addEnergyEvalLocal",__FILE__,__LINE__);
            }
            if (mp.setEnergy(energy,energyError)){
                mp.didEnergyEval();
            }
        } else {
            silosForKey(s).addEnergyEvalLocal(s,p,energy,energyError);
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
            silosForKey(owner).didLocalPublish(owner,newPoint,s,0);
            foreach(p,mp;localPoints){
                auto newPos=refPos.dynVars.dVarStruct.emptyX(true);
                newPos[]=pos.x;
                if (mp!is null){
                    mp.checkIfNeighbor(newPos,newPoint,pSize);
                }
                newPos.giveBack();
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.publishPoint(s,owner,newPoint,pos,pSize,flags);
            });
            silosForKey(owner).didLocalPublish(owner,newPoint,s,1);
        } else {
            silosForKey(s).publishPoint(s,owner,newPoint,pos,pSize,flags);
        }
    }
    void publishedLocalPoint(SKey s,Point p){
        if (hasKey(s)){
            notifyLocalObservers(delegate void(ExplorationObserverI!(T) obs){
                obs.publishedLocalPoint(_key,p);
            });
        } else {
            silosForKey(s).publishedLocalPoint(s,p);
        }
    }
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to silos
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin eAndMin){
        if (hasKey(s)){
            foreach (nP;neighbors){
                auto mp=localPoints[nP];
                mp.addEnergyEvalOther(eAndMin);
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.neighborHasEnergy(s,neighbors,eAndMin);
            });
        } else {
            silosForKey(s).neighborHasEnergy(s,neighbors,eAndMin);
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors might be restricted to silos or not
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin eAndMin){
        if (hasKey(s)){
            foreach (nP;neighbors){
                auto mp=localPoints[nP];
                mp.addGradEvalOther(p,eAndMin);
            }
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.neighborHasGradient(s,p,neighbors,eAndMin);
            });
        } else {
            silosForKey(s).neighborHasGradient(s,p,neighbors,eAndMin);
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
                    auto ppPtr=point in publishingPoints;
                    if (ppPtr!is null) pp=*ppPtr;
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
    void didLocalPublish(SKey s,Point p0,SKey source,int level){
        if (s==SKeyVal.All){
            foreach(sK,sL;silos){
                sL.didLocalPublish(sK,p0,source,level);
            }
        } else if (hasKey(s)){
            PointBcastProgress!(T) pp;
            synchronized(publishingPoints){
                pp=publishingPoints[p0];
            }
            assert(pp!is null);
            pp.communicatePublish(source,level);
            notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
                obs.didLocalPublish(s,p0,source,level);
            });
        } else {
            silosForKey(s).didLocalPublish(s,p0,source,level);
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
    PointEMin[] energyForPointsLocal(SKey s,Point[] points,PointEMin[] ens){
        if (hasKey(s)){
            if (ens.length<points.length)
                ens.length=points.length;
            foreach(i,p;points){
                auto mp=localPoints[p];
                ens[i]=mp.pointEMin;
            }
            return ens;
        } else {
            return silosForKey(s).energyForPointsLocal(s,points,ens);
        }
    }
    /// energy for the points (NAN if not yet known)
    PointEMin[] energyForPoints(SKey s,Point[]points,PointEMin[] ens){
        if (hasKey(s)){
            if (ens.length<points.length)
                ens.length=points.length;
            auto pO=new PointToOwner(&this.ownerOfPoint,delegate Point(size_t i){
                return points[i];
            }, points.length,false);
            Real[128] buf;
            auto lMem=LocalMem(buf);
            foreach (iter;pO){
                auto enLoc=lMem.allocArr!(PointEMin)(iter.localUb-iter.localLb);
                auto enLoc2=energyForPointsLocal(iter.owner,iter.points,enLoc);
                foreach(i,j;iter.idx){
                    ens[j]=enLoc2[i];
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
    MainPointI!(T) newPointAt(DynPVector!(T,XType) newPos,Point proposedPoint){
        MainPoint!(T) newP=pPool.getObj();
        refPos.checkX();
        auto p2=refPos.dup(PSDupLevel.EmptyDyn);
        p2.checkX();
        newP.pos.checkX();
        newP.pos.dynVars.x[]=newPos;
        synchronized(localPoints){ // don't allow null mainpoints...
            Point np;
            synchronized(owner){ // generate a point with a small propability of collisions
                bool setPoint=false;
                if (proposedPoint.isValid){
                    synchronized(localPointsKeys){
                        if (localPointsKeys.length==0 || proposedPoint.data>localPointsKeys[localPointsKeys.length-1].data){
                            localPointsKeys.appendEl(proposedPoint);
                            np=proposedPoint;
                            setPoint=true;
                        }
                    }
                }
                while(!setPoint){
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
        synchronized(this){
            insertAt(observers,ObserversEl(o));
        }
    }
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o){
        ObserversEl* newList;
        ObserversEl* newListLast;
        synchronized(this){
            auto pAtt=observers;
            while(pAtt!is null && pAtt.val!is o){
                newListLast=ObserversEl(pAtt.val,newList);
                if (newList is null) newList=newListLast;
            }
            if (pAtt.val is o){
                if (newListLast!is null) {
                    newListLast.next=pAtt.next;
                    writeBarrier();
                    observers=newList;
                } else {
                    writeBarrier();
                    observers=pAtt.next;
                }
            } else {
                assert(pAtt is null);
                while(newList!is null){
                    newListLast=newList.next;
                    delete newList;
                    newList=newListLast;
                }
            }
        }
    }
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o){
        synchronized(explorers){
            explorers[o.name]=o;
        }
        synchronized(generatingExplorers){
            generatingExplorers[o.name]=o;
        }
    }
    /// removes the given explorer
    bool rmExplorerNamed(char[] name){
        synchronized(generatingExplorers){
            return generatingExplorers.removeKey(name);
        }
    }
    void notifyLocalObservers(void delegate(ExplorationObserverI!(T))op){
        auto oAtt=observers;
        while(oAtt!is null){
            op(oAtt.val);
            oAtt=oAtt.next;
        }
    }
    /// adds a component
    void addComponent(SilosComponentI!(T)component){
        synchronized(_components){
            if (component.name in _components){
                if (_components[component.name]!is component){
                    throw new Exception("collision with component named "~component.name,__FILE__,__LINE__);
                }
            } else {
                _components[component.name]=component;
            }
        }
    }
    /// removes the component with the given name (stopping it)
    bool rmComponentNamed(string name){
        synchronized(_components){
            return _components.removeKey(name);
        }
    }
    /// returns the actual components
    HashMap!(string,SilosComponentI!(T)) components(){
        return _components;
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
        if (pC !is null && cachedP.lastSync>t){
            cachedP.mainPoint.retain;
            return cachedP.mainPoint;
        }
        cachedP.mainPoint=pPool.getObj();
        cachedP.mainPoint.pos.checkX();
        cachedP.mainPoint._point=p;
        cachedP.mainPoint._gFlags|=GFlags.LocalCopy;
        cachedP.lastSync=ev_time();
        cachedP.mainPoint[]=mainPoint(ownerOfPoint(p),p);
        synchronized(localCache){
            localCache[p]=cachedP;
        }
        cachedP.mainPoint.retain;
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
            auto lP=p.point in localCache;
            if (lP!is null && lP.mainPoint is p){
                p.release;
            }
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
        this._localPointsKeys=new BatchedGrowableArray!(Point,batchSize)();
        this.owner=new HashMap!(Point,SKey)();
        this.localCache=new HashMap!(Point,CachedPoint!(T))();
        this.publishingPoints=new HashMap!(Point,PointBcastProgress!(T))();
        this.pendingEvals=new PendingEvals!(T)(this);
        this.localOps=new Deque!(EvalOp!(T))();
        this.localPoints=new DenseLocalPointArray!(MainPointI!(T))(this.localPointsKeys); /// local points (position,...)
        this._rand=new RandomSync();
        this._refPos=refPos;
        this._constraints=constraints;
        this.log=log;
        this._nCenter=nCenter;
        if (nCenter is null) this._nCenter=new NotificationCenter();
        this.explorers=new HashMap!(char[],ExplorerI!(T))();
        this.generatingExplorers=new HashMap!(char[],ExplorerI!(T))();
        this.observers=null;
        this._components=new HashMap!(string,SilosComponentI!(T))();
        this.evaluator=cast(Method)this.input.evaluator.contentObj;
        if (this.input.evaluatorTask!is null) {
            this.evaluatorTask=cast(Sampler)this.input.evaluatorTask.contentObj;
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
                generatingExplorers[expl.name]=expl;
                addObserver(expl);
            } else {
                addObserver(observerT!(T)(cast(ExplorationObserverGen)o,this));
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
	logMsg(delegate void(CharSink s){
		dumper(s)("url: ")(silosCoreUrl());
	    });
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
    
    /// loops on all the local points
    int opApply(int delegate(ref Point,ref MainPointI!(T) el)loopBody){
        return localPoints.opApply(loopBody);
    }
    
}

private PNetSilos!(Real) dummyPNetSilos;
