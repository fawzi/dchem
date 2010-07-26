module dchem.pnet.PNetSilos;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import dchem.input.RootInput;
import dchem.pnet.PNetModels;
import blip.core.Variant;

struct LoadStats{
    Real load;
    EKey explorer;
    mixin(serializeSome("LoadStats","load|explorer"));
    mixin printOut!();
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
    Silos!(Real)[] silosReal;
    Silos!(LowP)[] silosLowP;
    
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

    // 
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
    Real dirSize=1.0;
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
    Real minRealNormSelf2{
        return cartesianDiffMin/2;
    }
    /// explorationSize
    Real explorationSize(){
        discretizationStep;
    }
    /// square of maximum distance from the optimal direction point in cartesian units
    Real dirCartesianSize2(){
        return cartesianDiffMin*cartesianDiffMin;
    }
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.Silos",`
        evaluator: the method to perform energy evaluations
        precision: the precision with which the points are stored
        explorers: the observers and explorers active on the point network
        loaders: SilosWorkers that load already evaluated points (called after the explorer are set up, before exploring)
        monitors: SilosWorkers called at regular intervals during the exploration
        finisher: list of finishers, SilosWorkers called after the exploration is finished
        minProjectionResidual: minimum residual after projecting out the linear dependend directions (0.05)
        sameDirCosAngle: minimum cos of then angle in dual space to consider a direction equivalent (0.866)
        minNormDual: minimum norm in the dual T space to accept an exploration (in units of explorationSize) (0.4)
        minNormDualSelf: minimum norm in the dual T space to accept an exploration for a self generated direction (in units of discretizationStep) (0.1)
        discretizationStep: discretization step (0.2)
        cartesianDiffMin: minimum cartesian difference (0.02)
        maxNormDual: max norm in the dual T space to accept an exploration (in units of explorationSize) (1.9)
        dirSize: square of the maximum distance from the optimal direction point in explorationSize units (1.0)
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
            auto silos=new Silos!(LowP)(this);
            synchronized(this){
                silosLowP~=silos;
            }
            silos.run();
        } else if (precision=="Real"){
            auto silos=new Silos!(Real)(this);
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
        if (c in splitChars){
            if ((!skipEmpty)||j>i){
                res~=str[i..j];
            }
            i=j+1;
        }
    }
    if (i<str.length) res~=str[i..$];
    return res;
}

char[] realFromInput(char[]props){
    char[] res;
    foreach(property;ctfeSplit("| \n",props,true){
        if (property.length>0){
            res~=`
            T `~property~`(){
                return cast(T)input.`~property~`;
            }`;
        }
    }
    return res;
}

class Silos(T): ExplorerI!(T){
    SilosGen input;
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    EKey key;
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
    /// reference position/particle system
    ParticleSys!(T) _refPos;
    /// constraints (taken from the evaluator, guaranteed to be non null)
    ConstraintI!(T) _constraints;
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
    
    mixin(realFromInput(`discretizationStep|minProjectionResidual|sameDirCosAngle|minNormDual|minNormDualSelf
        minRealNormSelf0|minRealNormSelf1|maxNormDual|explorationSize|dirSize2|dirCartesianSize2`));
    
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
        cInstance=getInstanceForClass(InstanceGetFlags.ReuseCache|InstanceGetFlags.NoAllocSubOpt|InstanceGetFlags.Wait);
        switch(cInstance.activePrecision()){
        case Precision.Real:
            _refPos=cInstance.pSysReal.dupT!(T)();
            _constraints=cInstance.constraintsReal();
            break;
        case Precision.LowP:
            _refPos=cInstance.pSysLowP.dupT!(T)();
            _constraints=cInstance.constraintsLowP();
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
    
    MainPoint hydrate(LocalSilosI!(T)localContext){
        auto pSys=localContext.refPos.dup();
        pSys[]=pos;
        DynPVect!(T,2) mDir;
        if (minDir.isNonNull()){
            mDir=minDir.toPVector(pSys.dualDxGroup,true);
        }
        MainPoint res=new MainPoint(localContext,point,pSys,gFlags|GFlags.LocalCopy,explorationSize,
            mDir,minDirScale,exploredDirs,neighbors);
        return res;
    }
    
}

