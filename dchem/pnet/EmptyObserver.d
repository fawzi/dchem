module dchem.pnet.EmptyObserver;
import dchem.Common;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;
import dchem.sys.ParticleSys;
import dchem.input.WriteOut;

/// empty observer (just as superclass to simplify definition of new observers)
class EmptyObserver(T): ExplorationObserverI!(T) {
    this(){}
    /// increases the runLevel on all silos, i.e. you should call it only with SKeyVal.All
    /// (at the moment there is no support for dynamic adding/removal of silos)
    void increaseRunLevel(SKey s,RunLevel speed){}
    /// adds energy for a point local to s and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy,Real energyError){}
    /// adds gradient value to a point that should be owned by s. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){}
    /// communicates to s that the given point is being explored
    /// pSize is the point size, flags the flags of the point
    void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){}
    /// communicates that the given local point has been successfully published
    void publishedLocalPoint(SKey s,Point point){}
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to s
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin energy){}
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin energy){}
    
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point,SKey owner){}
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source,int level){}
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey,Point){}
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradientLocal(SKey s,Point p,Real energy,Real energyError){
        return false;
    }
    /// checks it local point is somehow invalid and should better be skipped
    bool shouldFilterLocalPoint(SKey s,Point p){
        return false;
    }
    /// unique name to identify this observer (the different processes should use the same name)
    char[] name(){
        assert(0,"to be implemented in the subclasses"); // implement??
    }
    
}
