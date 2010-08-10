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

struct LoadStats{
    Real load;
    SKey silosKey;
    mixin(serializeSome("LoadStats","load|silosKey"));
    mixin printOut!();
    int opCmp(LoadStats l2){
        return ((load<l2.load)?-1:((load==l2.load)?cmp(cast(long)silosKey,cast(long)l2.silosKey):1));
    }
}
/+
class RemotePointEval(T):RemoteTask{
    Point point;
    char[] ownerUrl;
    CalculationContext ctx;
    mixin(serializeSome("dchem.minEE.RemotePointEval!("~T.stringof~")","point|ownerUrl"));
    
    this(){}
    this(Point p,char[]oUrl){
        point=p;
        ownerUrl=oUrl;
    }
    void execute(Variant args){
        auto ctx=cast(ActiveExplorer!(T))cast(Object)ProtocolHandler.proxyForUrl(ownerUrl);
        auto pPos=ctx.pointPos(point); // get real point
        ctx=args.get!(CalculationContext)();
    }
    void stop(){
        if (ctx!is null) ctx.stop();
        ctx=null;
    }
    
}
class MasterCalculator:Sampler{
    InputField method;
    mixin myFieldMixin!();

    ParticleSys!(Real) refPointReal;
    ParticleSys!(LowP) refPointLowP;

    Set!(PointAndDir) inExploration; /// point directions that are in exploration
    Set!(Point) inEvaluation; /// points that are in evaluation
    Set!(Point) toEvaluate; /// points that should be avaluated
    size_t overPrepare=0; /// number of extra calculation to prepare
    ExplorerI!(T)[] explorers; /// explorers that can create new work
    
    /// adds a context that can calculate something
    void addContext(CalculationContext c){
        synchronized(this){
            waitingContexts.append(c);
        }
        update();
    }
    /// adds a point to evaluate coming from PointAndDir
    void addPointForPointDir(Point p,PointAndDir origin){
        synchronized(this){
            if (p.isValid){
                toEvaluate.add(p);
            }
            inExploration.remove(origin);
        }
        update();
    }
    /// updates the calculator: tries to give work to waitingContexts and request more work if needed
    void update(){
        Point p;
        CalculationContext ctx;
        while(1){
            ctx=null;
            synchronized(this){
                if (!toEvaluate.take(p)) break;
                if (!waitingContexts.popFront(ctx)) break;
            }
            // ctx.remoteExe...
        }
        if (p.isValid && ctx is null) {
            synchronized(this){
                toEvaluate.add(p);
            }
        }
    }
}+/

class SilosGen:Sampler{
    PNetSilos!(Real)[] silosReal;
    PNetSilos!(LowP)[] silosLowP;
    
    /// number of exploration steps to perform
    long explorationSteps;
    /// number of exploration steps to perform
    long maxExplorationTries=10;
    /// the current method for energy evaluation
    InputField evaluator;

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
        if (precision!="LowP" && precision!="Real"){
            log("precision should be either LowP or Real, not '")(precision)("' in field ")(myField)("\n");
            res=false;
        }
        foreach(e;explorers){
            auto o=cast(Object)e.content;
            if (cast(ExplorationObserverGen)o is null){
                log("explorers should be of either explorers or observer, not '")(o.classinfo.toString())("' in field ")(myField)("\n");
                res=false;
            }
        }
        foreach(l;loaders){
            auto o=cast(Object)l.content;
            if (cast(SilosWorkerGen)o is null){
                log("loaders should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myField)("\n");
                res=false;
            }
        }
        foreach(m;monitors){
            auto o=cast(Object)m.content;
            if (cast(SilosWorkerGen)o is null){
                log("monitors should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myField)("\n");
                res=false;
            }
        }
        foreach(f;finishers){
            auto o=cast(Object)f.content;
            if (cast(SilosWorkerGen)o is null){
                log("finishers should be a SilosWorker not '")(o.classinfo.toString())("' in field ")(myField)("\n");
                res=false;
            }
        }
        return res;
    }

    void run(){
        if (precision=="LowP"){
            auto silos=new PNetSilos!(LowP)(this);
            synchronized(this){
                silosLowP~=silos;
            }
            silos.run();
        } else if (precision=="Real"){
            auto silos=new PNetSilos!(Real)(this);
            synchronized(this){
                silosReal~=silos;
            }
            silos.run();
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

char[][] ctfeSplit(char[] splitChars,char[]str,bool skipEmpty){
    char[][]res;
    size_t i=0;
    foreach(j,c;str){
        foreach (c2;splitChars){
            if (c==c2){
                if ((!skipEmpty)||j>i){
                    res~=str[i..j];
                }
                i=j+1;
            }
        }
    }
    if (i<str.length) res~=str[i..$];
    return res;
}

char[] realFromInput(char[]props){
    char[] res;
    foreach(property;ctfeSplit("| \n",props,true)){
        if (property.length>0){
            res~=`
            T `~property~`(){
                return cast(T)input.`~property~`;
            }`;
        }
    }
    return res;
}

struct CachedPoint(T){
    MainPoint!(T) mainPoint;
    Time lastSync;
}

class PNetSilos(T): LocalSilosI!(T){
    enum RunLevel:int{
        Setup,
        Running,
        Stop
    }
    RunLevel runLevel;
    SilosGen input;
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    SKey _key;
    SKey key(){ return _key; }
    HashMap!(SKey,PNetSilosI!(T)) silos;
    MinHeapSync!(LoadStats) loads; 
    BatchedGrowableArray!(Point,batchSize) localPointsKeys; // indexes for local owned points
    HashMap!(Point,SKey) owner;
    HashMap!(Point,CachedPoint!(T)) localCache;
    CalculationContext[Point] localCalcInProgress;
    Set!(Point) calcInProgress;
    DenseLocalPointArray!(MainPointI!(T)) localPoints; /// local points (position,...)
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
    Deque!(ExplorerI!(T)) explorers;
    /// observers
    Deque!(ExplorationObserverI!(T)) observers;
    /// evaluator
    Method evaluator;
    PoolI!(MainPoint!(T)) pPool;
    
    MainPoint!(T)allocPoint(PoolI!(MainPoint!(T))p){
        return new MainPoint!(T)(this,p);
    }
    this(SilosGen sGen){
        input=sGen;
        pPool=cachedPool(&this.allocPoint);
        assert(0,"to do: complete init");
    }
    /// adds the computed gradient to the point p
    void addGradEval(Point p,PSysWriter!(T) pSys){
        auto mp=localPoints[p];
        synchronized(mp){
            auto e=mp.pos.dynVars.potentialEnergy;
            mp.pos[]=pSys;
            if (isNaN(mp.pos.dynVars.potentialEnergy)) mp.pos.dynVars.potentialEnergy=e;
        }
        mp.didGradEval();
    }
    /// writes out a log message
    void logMsg(void delegate(CharSink)writer){
        char[512] buf;
        auto gArr=lGrowableArray(buf,0);
        auto s=dumper(&gArr.appendArr);
        s("<MinEELog id=\"")(key)("\" time=\"")(Clock.now.ticks)("\">");
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
    bool speculativeGradient(Point p,Real energy){
        bool res=false;
        int i=0;
        foreach(e;explorers){
            res=res || e.speculativeGradient(p,energy);
        }
        struct EAdd{
            void delegate(Point,Real) op;
            Point p;
            Real energy;
            void doOp(){
                op(p,energy);
            }
        }
        auto eAdd=new EAdd;
        eAdd.energy=energy;
        eAdd.p=p;
        eAdd.op=&addEnergyEvalLocal;
        Task("addOp",&eAdd.doOp).autorelease.submit(defaultTask);
        return res;
    }
    
    mixin(realFromInput(`discretizationStep|minProjectionResidual|sameDirCosAngle|minNormDual|minNormDualSelf
        minRealNormSelf0|minRealNormSelf1|maxNormDual|explorationSize|dirSize2|dirCartesianSize2
        maxMoveInternal|maxMoveCartesian|maxMoveDual|maxNormCartesian|minMoveInternal|sameDirCosAngleCartesian|
        minRealNormSelf2|dirDualSize2|inDirCartesianScale2|zeroLen`));
    
    /// returns a globally unique string 
    char[] nextUniqueStr(){
        return collectAppender(delegate void(CharSink s){
            s(nameId); s("_"); writeOut(s,nextLocalId.next());
        });
    }
    /// returns a most likely valid point id
    Point nextPointId(){
        ushort r;
        rand()(r);
        return Point(((nextPntNr.next())<<12)|cast(ulong)(r&0xFFF));
    }
    
    /// returns true if the evaluation of the given point is in progress. Works only master/owner...?
    bool isInProgress(Point p){
        bool res=false;
        synchronized(this){
            res=calcInProgress.contains(p);
        }
        return res;
    }
    
    mixin(serializeSome("dchem.MinEExplorer_"~T.stringof,``));
    mixin printOut!();
    
    void run(){
        CalculationContext cInstance=evaluator.getCalculator(true,[]);
        switch(cInstance.activePrecision()){
        case Precision.Real:
            _refPos=cInstance.pSysReal.dupT!(T)(PSDupLevel.All);
            auto c1=cInstance.constraintsReal();
            if (c1 !is null)
                _constraints=constraintT!(T)(c1.constraintGen,_refPos);
            break;
        case Precision.LowP:
            _refPos=cInstance.pSysLowP.dupT!(T)(PSDupLevel.All);
            auto c2=cInstance.constraintsLowP();
            if (c2 !is null)
                _constraints=constraintT!(T)(c2.constraintGen,_refPos);
            break;
        default:
            assert(0,"unknown activePrecision");
        }
        if (_constraints is null) _constraints=new NoConstraint!(T)();
        // run the loaders
        //loaders
        // if empty & no explorator, add point at actual evaluator point
        /+
        while(1)
            get context
            choose eval host
            remote task on it

        remote task:
            start independent task that executes it
            while 1
                getNext point from master if invalid, brek
                decide if eval
                if eval
                    eval, return res,break

        +/
    }
    
    bool nextExplorer(ref ExplorerI!(T) res){
        synchronized(explorers){
            if (explorers.popFront(res)){
                explorers.append(res);
            } else {
                res=null;
            }
        }
        return res!is null;
    }
    /// returns the next point to evaluate
    Point pointToEvaluate(){
        ExplorerI!(T) e;
        for (int i=0;i<input.maxExplorationTries;++i){
            if (!nextExplorer(e)){
                break;
            }
            auto res=e.pointToEvaluate();
            if (res.isValid) return res;
        }
        return Point(0); // no valid next computation
    }
    
    /// returns the key of the silos with the lowest expected load
    SKey nextFreeSilos(){
        auto w=loads.popWait();
        mixin(mkActionMixin("updLoad","w|silos|loads",`
        w.load=silos[w.silosKey].load();
        loads.push(w);`));
        if (key==w.silosKey){ // local update
            updLoad();
        } else {
            Task("updLoad",updLoad).autorelease.submitYield();
        }
        return w.silosKey;
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
    void addEnergyEvalLocal(Point p,Real energy){
        auto mp=localPoints[p];
        if (! mp.isLocalCopy){
            throw new Exception("Local copy in PNetSilos.addEnergyEvalLocal",__FILE__,__LINE__);
        }
        mp.addEnergyEvalLocal(energy);
        foreach(obs;observers){
            obs.addEnergyEvalLocal(p,energy);
        }
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(Point p,PSysWriter!(T) pSys){
        auto mp=localPoints[p];
        if (! mp.isLocalCopy){
            throw new Exception("Local copy in PNetSilos.addGradEvalLocal",__FILE__,__LINE__);
        }
        mp.pos[]=pSys;
        mp.didGradEval();
        foreach(obs;observers){
            obs.addGradEvalLocal(p,pSys);
        }
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void addExploredPoint(SKey owner,Point newPoint,PSysWriter!(T) pos,T pSize,uint flags){
        this.owner[newPoint]=owner;
        foreach(p,mp;localPoints){
            auto newPos=refPos.dynVars.dVarStruct.emptyX();
            newPos[]=pos.x;
            if (mp!is null){
                mp.checkIfNeighbor(newPos,newPoint,pSize);
            }
            newPos.giveBack();
        }
        foreach(obs;observers){
            obs.addExploredPoint(owner,newPoint,pos,pSize,flags);
        }
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to silos
    void neighborHasEnergy(Point p,Point[] neighbors,Real energy){
        foreach (nP;neighbors){
            auto mp=localPoints[nP];
            mp.addEnergyEvalOther(p,energy);
        }
        foreach(obs;observers){
            obs.neighborHasEnergy(p,neighbors,energy);
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors might be restricted to silos or not
    void neighborHasGradient(LazyMPLoader!(T)p, Point[] neighbors, Real energy){
        foreach (nP;neighbors){
            auto mp=localPoints[nP];
            mp.addGradEvalOther(p,energy);
        }
        foreach(obs;observers){
            obs.neighborHasGradient(p,neighbors,energy);
        }
    }
    
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(LocalSilosI!(T) silos,Point point){
        foreach(obs;observers){
            obs.finishedExploringPointLocal(silos,point);
        }
    }
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(Point point){
        foreach(obs;observers){
            obs.finishedExploringPoint(point);
        }
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void dropPoint(SKey sKey,Point point){
        if (sKey==key){
            auto mp=localPoints[point];
            mp.drop();
        }
        foreach(obs;observers){
            obs.dropPoint(sKey,point);
        }
    }
    
    // PNetSilosI
    
    /// stops all silos (at the moment there is no support for dynamic adding/removal of silos, worth adding???)
    void shutdown(){
        bool shouldStop=false;
        synchronized(this){
            if (runLevel<RunLevel.Stop){
                runLevel=RunLevel.Stop;
                shouldStop=true;
            }
        }
        if (shouldStop){ // square b-cast, change??
            foreach(k,sil;silos){
                sil.shutdown();
            }
        }
    }
    /// load (usage) of the silos in some arbitrary units
    Real load(){
        return 0.001*localPoints.length;
    }
    /// returns the position of the given local point (remove??)
    PSysWriter!(T)pointPosLocal(Point p){
        auto mp=localPoints[p];
        return pSysWriter(mp.pos);
    }
    /// returns the position of the given point (remove??)
    PSysWriter!(T)pointPos(Point p){
        auto ownr=ownerOfPoint(p);
        if (ownr==key){
            return pointPosLocal(p);
        }
        return silosForKey(ownr).pointPosLocal(p);
    }
    /// energy for the local points (NAN if not yet known)
    Real[] energyForPointsLocal(Point[] points,Real[] ens){
        if (ens.length<points.length)
            ens.length=points.length;
        foreach(i,p;points){
            auto mp=localPoints[p];
            ens[i]=mp.energy;
        }
        return ens;
    }
    /// energy for the points (NAN if not yet known)
    Real[] energyForPoints(Point[]points,Real[] ens){
        if (ens.length<points.length)
            ens.length=points.length;
        auto pO=new PointToOwner(&this.ownerOfPoint,delegate Point(size_t i){
            return points[i];
        }, points.length,false);
        Real[128] buf;
        auto lMem=LocalMem(buf);
        foreach (iter;pO){
            auto kAtt=iter.owner;
            if (kAtt==key){
                auto enLoc=lMem.allocArr!(Real)(iter.localUb-iter.localLb);
                enLoc=energyForPointsLocal(iter.points,enLoc);
                foreach(i,j;iter.idx){
                    ens[j]=enLoc[i];
                }
                lMem.deallocArr(enLoc);
            }
        }
        return ens;
    }
    /// returns a snapshot of the given point
    DriedPoint!(T)mainPoint(Point p){
        auto ownr=ownerOfPoint(p);
        if (ownr==key){
            return mainPointLocal(p);
        }
        return silosForKey(ownr).mainPointLocal(p);
    }
    /// returns a snapshot of the given point that is local to this silos
    DriedPoint!(T)mainPointLocal(Point p){
        auto mp=localPoints[p];
        return mp.driedPoint;
    }
    /// tells the local points neighs that they have p0 as neighbor
    void addPointToLocalNeighs(Point p0,Point[]neighs){
        foreach(np;neighs){
            auto mp=localPoints[np];
            mp.addNeighbor(p0);
        }
    }
    /// tells the local point p0 that neighDirs (computed when p0 gradient was hadGrad) should be added to it
    bool addNeighDirsToLocalPoint(Point p0,PointAndDir[]neighDirs,DirDistances!(T)[]dirDists,bool hadGrad){
        auto mp=localPoints[p0];
        return mp.addNeighbors(neighDirs,dirDists,hadGrad);
    }
    /// informs that source has processed point p0
    void processedLocal(Point p0,SKey source){
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
        explorers.append(o);
    }
    /// removes the given explorer
    void rmExplorer(ExplorerI!(T)o){
        explorers.filterInPlace(delegate bool(ExplorerI!(T) i){ return i!is o; });
    }
    SKey ownerOfPoint(Point p){
        return owner[p];
    }
    PNetSilosI!(T) silosForKey(SKey k){
        if (k==key) return this;
        return silos[k];
    }
    /// local point mainpoint (the real reference point)
    MainPointI!(T) mainPointL(Point p){
        return localPoints[p];
    }
    /// creates a new point in this silos located at newPos, the point is not yet broadcasted
    Point newLocalPointAt(DynPVector!(T,XType) newPos){
        auto newPoint=nextPointId();
        auto newP=pPool.getObj();
        newP.pos.dynVars.checkX();
        newP.pos.dynVars.x[]=newPos;
        localPoints[newPoint]=newP;
        return newPoint;
    }
    /// a local point (possibly a copy), is retained, and needs to be released (thus the create in the name)
    /// the time t is used to decide if a cached version can be used
    MainPointI!(T)createLocalPoint(Point p,Time t){
        if (ownerOfPoint(p)==key){
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
        cachedP.lastSync=Clock.now;
        cachedP.mainPoint[]=mainPoint(p);
        synchronized(localCache){
            localCache[p]=cachedP;
        }
        return cachedP.mainPoint;
    }
    /// makes a point "public" informing other silos that that region has been explored
    void bcastPoint(Point p){
        auto mp=localPoints[p];
        foreach(s;silos){
            s.addExploredPoint(key,p,pSysWriter(mp.pos),mp.explorationSize,mp.gFlags);
        }
    }
    /// drops a cached point (the point is not in use anymore)
    void dropCachedPoint(MainPointI!(T)p){
        synchronized(localCache){
            localCache.removeKey(p.point);
        }
    }
    
}

