module dchem.pnet.MinEExplorer;
import dchem.pnet.PNetModels;

/// a structure to keep point and its energy together
struct PointAndEnergy{
    Point point;
    Real energy;
    mixin(serializeSome("PointAndEnergy","point|energy"));
    mixin printOut!();
}

class MinEExplorerDef:InputElement{
    long nEval;
    char[] precision;
    mixin myFieldMixin!();
    mixin(SerializeSome("dchem.MinEExplorer",`
    nEval : number of evaluations to perform`));
    
    bool verify(CharSink s){
        return true;
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
class MinEExplorer(T):ExplorerI(T){
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    HashSet!(Point) removedPoints;
    // ExplorationObserverI(T)
    
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(LocalSilosI!(T) silos,Point p,Real energy){
        toExploreMore.push(PointAndEnergy(p,energy));
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(LocalSilosI!(T) silos,Point p,PSysWriter!(T) pSys){
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void addExploredPoint(SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
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
    void dropPoint(SKey e,Point p){
        rmPoint(p);
    }
    
    // ExplorerI(T)
    
    /// returns a point to evaluate
    Point pointToEvaluate(){
        while {
            auto pe=toExploreMore.pop();
            if (!removedPoints.contains(pe.point)){
                return pe.point;
                /+pDir=lowestPoint.exploreNext;
                find silos
                pOld=silos.localPoint(pDir.point);
                pNew=pOld.createPosInDir(pDir);
                newP=silos.newPoint(pNew);
                if(pOld.acceptableDir(newP)){
                	bCast(newP);
                	return newP to evaluate;
                } else {
                	drop(newP);
                }
                try again +/
            }
        }
    }
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(Point,uint flags){
        rmPoint(p);
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradient(Point p,Real energy){
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
}
