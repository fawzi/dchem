module dchem.pnet.MinEExplorer;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;

class MinEExplorerDef:SilosWorkerGen{
    long nEval;
    char[] precision;
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.MinEExplorer",`
    nEval : number of evaluations to perform`));
    
    bool verify(CharSink s){
        return true;
    }
    
    SilosWorkerI!(Real) silosWorkerReal(){
        auto res=new MinEExplorer!(Real)(this);
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        auto res=new MinEExplorer!(LowP)(this);
    }
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

/// an object that can offer new points to explore
/// actually should not inherit from ExplorationObserverI, but this way we avoid multiple inheritance bugs
class MinEExplorer(T):ExplorerI!(T),SilosWorkerI!(T){
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    HashSet!(Point) removedPoints;
    void delegate(ExplorerI!(T)) available;
    MinEExplorerDef input;
    this(MinEExplorerDef input){
        this.input=input;
        this.toExploreMore=new MinHeapSync!(PointAndEnergy)();
        this.removedPoints=new HashSet!(Point);
    }
    
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(LocalSilosI!(T) silos,Point p,Real energy){
        synchronized(this){
            toExploreMore.push(PointAndEnergy(p,energy));
            if (available!is null){ // protect against false calls?
                available(this);
                available=null;
            }
        }
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(LocalSilosI!(T) silos,Point p,PSysWriter!(T) pSys){
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void publishPoint(SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
        if (!isNAN(pos.potentialEnergy)){
            assert((flags&MainPoint!(T).GFlags.EnergyInfo)=MainPoint!(T).GFlags.EnergyKnown,"non NAN energy, but not EnergyKnown");
            if ((flags&(MainPoint!(T).GFlags.DoNotExplore|MainPoint!(T).GFlags.FullyExplored))==0){
                toExploreMore.add(p);
            }
        }
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    void neighborHasEnergy(Point p,Point[] neighbors,Real energy){ }
    /// a neighbor point has calculated its gradient (and energy)
    void neighborHasGradient(Point p,Point[] neighbors, PSysWriter!(T)pSys){ }
    
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(LocalSilosI!(T) silos,Point p){
        rmPoint(p);
    }
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(Point p){
        rmPoint(p);
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey e,Point p){
        rmPoint(p);
    }
    
    // ExplorerI(T)
    
    /// returns a point to evaluate
    Point pointToEvaluateLocal(SKey k,void delegate(ExplorerI!(T))available){
        bool availableCalled=false;
        while(true) {
            PointAndEnergy pe;
            if (toExploreMore.popNext(pe)){
                if (!availableCalled) {
                    available();
                    availableCalled=true;
                }
                if (!removedPoints.contains(pe.point)){
                    auto lowestPoint=silos.createLocalPoint(pe);
                    auto pDir=lowestPoint.exploreNext();
                    if (pDir.point.data<=1){
                        continue;
                    }
                    toExploreMore.add(pe);
                    auto pOld=silos.localPoint(pDir.point);
                    auto posNew=pOld.createPosInDir(pDir);
                    newP=silos.newPoint(posNew);
                    if(pOld.acceptableDir(newP)){
                        bCast(newP);
                        return newP;
                    } else {
                        drop(newP);
                    }
                    //try again
                }
            } else {
                assert(this.available==null || this.available is available,"invalid available value");
                this.available=available;
                return Point(1); // needs to wait...
            }
        }
    }
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(Point,uint flags){
        rmPoint(p);
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradient(SKey s,Point p,Real energy){
        synchronized(toExploreMore){
            auto l=toExploreMore.length;
            if (l>0 && heap.data[0].energy>=energy){
                return true;
            }
        }
        return false;
    }
    void rmPoint(Point p){
        removedPoints.add(p);
    }
    
    void workOn(LocalSilosI!(T) silos){
        silos.addExplorer(this);
    }
    
}
