module dchem.pnet.PNetSilos;


struct LoadStats{
    Real load;
    EKey explorer;
    mixin(serializeSome("LoadStats","load|explorer"));
    mixin printOut!();
}

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
class MasterCalculator{
    Set!(PointAndDir) inExploration; /// point directions that are in exploration
    Set!(Point) inEvaluation; /// points that are in evaluation
    Set!(Point) toEvaluate; /// points that should be avaluated
    Deque!(CalculationContext) waitingContexts; /// contexts that can be used to do evaluations
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
}

class MinEExplorer(T): Sampler,ActiveExplorer!(T){
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    EKey key;
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    /// list of all the currently active explorers (including this one)
    ActiveExplorer!(T)[EKey] activeExplorers;
    MinHeapSync!(LoadStats) loads; 
    BatchGrowableArray!(Point,batchSize) localPointsKeys; // indexes for local owned points
    HashMap!(Point,EKey) owner;
    CalculationInstance[Point] localCalcInProgress;
    Set!(Point) calcInProgress;
    DenseLocalPointArray!(MainPoint!(T)) localPoints; /// local points (position,...)
    Random randLocal;
    RandomSync rand;
    /// step used for the discretization
    T discretizationStep;
    /// place to store the trajectory (journal)
    char[] trajDir;
    /// the current method for energy evaluation
    InputField evaluator;
    /// constraints (taken from the evaluator)
    MultiConstraint _constraints;
    /// place to log messages
    CharSink log;
    /// notification center
    NotificationCenter nCenter;
    /// explorers
    ExplorerI!(T)[] explorers;
    /// observers
    ExplorationObserverI!(T)[] observers;
    
    /// adds the computed gradient to the point p
    void addGradEval(Point p,PSysWriter!(T) pSys){
        auto mp=localPoints[p];
        synchronized(mp){
            auto e=mp.pos.dynVars.potentialEnergy;
            mp.pos[]=pSys;
            if (isNAN(mp.pos.dynVars.potentialEnergy)) mp.pos.dynVars.potentialEnergy=e;
        }
        mp.didGradEval();
    }
    /// writes out a log message
    void logMsg(void delegate(CharSink)writer){
        char[512] buf;
        auto gArr=lGrowableArray(buf,0);
        auto s=dumper(&gArr.appendArr);
        s("<MinEELog id=\"")(key)("\" time=\"")(toString(Clock.now()))("\">");
        writer(s.call);
        s("</MinEELog>\n");
    }
    /// writes out a log message
    void logMsg(char[]msg){
        logMsg(delegate void(CharSink s){ s(msg); });
    }
    /// constraints of the current system
    MultiConstraint constraints(){
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
        while(1){
            e;explorers){
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
    
    /// minimum norm of uint vector after projection of null direction to be considered as a valid direction
    T minProjectionResidual(){
        return 0.01;
    }
    /// minimum cos of then angle in dual space to consider a direction equivalent
    T sameDirCosAngle(){
        return 0.866;
    }
    /// minimum norm in the dual T space to accept an exploration (in units of explorationSize)
    T minNormDual(){
        return 0.4;
    }
    /// minimum norm in the dual T space to accept an exploration for a self generated direction 
    /// (in units of explorationSize)
    T minNormDualSelf(){
        return 0.1;
    }
    /// minimum norm in the real (cartesian) space for a self generated direction to be accepted before moving
    T minRealNormSelf0(){
        return cartesianDiffMin/4;
    }
    /// minimum norm in the real (cartesian) space to which a self generated direction should be rescaled
    T minRealNormSelf1(){
        return cartesianDiffMin;
    }
    /// minimum norm in the real (cartesian) space for the real movement (after constraints,...)
    /// of a self generated direction to be accepted
    T minRealNormSelf2(){
        return cartesianDiffMin/2;
    }
    /// max norm in the dual T space to accept an exploration (in units of explorationSize)
    T maxNormDual(){
        return 1.9;
    }
    /// explorationSize
    T explorationSize(){ return discretizationStep; }
    /// square of the maximum distance from the optimal direction point in explorationSize units
    T dirSize2(){
        return explorationSize();
    }
    /// square of maximum distance from the optimal direction point in cartesian units
    T dirCartesianSize2(){
        return cartesianDiffMin();
    }
    
    /// repulsionSize
    T repulsionSize(){ return discretizationStep; }
    bool pointIsExplorable(EKey eK,Point p){
        MainPoint!(T) mainP;
        synchronized(this){ // sync could be avoided most of the time 
            mainP=localPoints[p];
        }
        if (mainP is null) throw new Exception("asked unknown point "~p.toString(),__FILE__,__LINE__);
        auto flags=mainP.gFlags;
        return (flags&GFlags.Evaluated)&&(!(flags&(GFlags.DoNotExplore|GFlags.FullyExplored|GFlags.FullyEvaluated|GFlags.OldApprox)));
    }
    /// returns a globally unique string 
    char[] nextUniqueStr(){
        return collectAppender(delegate void(CharSinker s){
            s(nameId); s("_"); writeOut(s,nextUniqueId.next());
        });
    }
    /// returns a most likely valid point id
    Point nextPointId(){
        ushort r;
        rand(r);
        return Point(((nextPntNr.next())<<12)|cast(ulong)(r&0xFFF));
    }
    
    /// returns true if the evaluation of the given point is in progress
    bool isInProgress(Point p){
        bool res=false;
        synchronized(this){
            res=(p in calcInProgress)!is null;
        }
        return res;
    }
    
    mixin(serializeSome("dchem.MinEExplorer_"~T.stringof,
        `trajDir: directory where to store the trajectory (journal)`));
    mixin printOut!();
    
    void run(){
        bool restarted=false;
        evaluator.method.setupCalculatorClass();
        // possiby restarts
        if (! restarted){
            cInstance=getInstanceForClass(InstanceGetFlags.ReuseCache|InstanceGetFlags.NoAllocSubOpt|InstanceGetFlags.Wait);
        }
        master.nextCalculation();
    }
    
    // exploration:
    // find smallest energy still "free"
    // evaluate next direction of it
    // apply constraints, if movement is too small declare it as fully visited and neighbor 1
    // bcast point as explored
    // possibly wait for non collision confirmation
    // start calculation
    // if too close to existing points stop calculation???
    // when calculation is finished
    // if mainpoint:
    //    gradient -> orient neighbors, compile visited flags, perform topology analysis (attractor, min,max,...)
    //    second deriv check
    // else store energy, first deriv check??
    
    /// evaluates the next computation
    PointAndDir nextComputation(){
        auto smallPoint=toExploreMore.waitPop();
        auto nextPDir=activeExplorers[owner[smallPoint.point]].exploreNext(cheapGrad);
        if (nextPDir.isValid && nextPDir.dir!=0){
            toExploreMore.push(smallPoint);
        } // else other explorers have a point that is not executable as first, but that should create no problems
        return nextPDir;
    }
    
    /// returns the key of the worker with the lowest expected load
    EKey findNextWorker(){
        auto w=loads.popWait();
        mixin(mkActionMixin("updLoad","w|activeExplorers|loads",`
        w.load=activeExplorers[w.key].load();
        loads.push(w);`));
        if (key==w.key){ // local update
            updLoad();
        } else {
            Task("updLoad",&updLoad).autorelease.submitYield();
        }
        return w.key;
    }
    /// evaluates the given point and dir, and either request a real evaluation there, or communicates that it is blocked
    Point evaluatePointAndDir(){
        activeExplorers[findNextWorker()];
    }
    
    void stop(){
        
    }
    
    bool verify(CharSink s){ return true; }
}
