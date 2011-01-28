module dchem.pnet.EmptyExplorer;
import dchem.Common;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;
import dchem.sys.ParticleSys;
import dchem.input.WriteOut;

/// an explorer that does nothing
/// not really useful, mainly for completness
class EmptyExplorerDef:SilosWorkerGen{
    this(){
    }
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.EmptyExplorer",``));
    mixin printOut!();
    bool verify(CharSink s){
        return true;
    }
    
    SilosWorkerI!(Real) silosWorkerReal(){
        auto res=new EmptyExplorer!(Real)();
        return res;
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        auto res=new EmptyExplorer!(LowP)();
        return res;
    }
}

/// an explorer that does nothing (useful as base class)
class EmptyExplorer(T):ExplorerI!(T),SilosWorkerI!(T){
    
    this(){
    }
    /// increases the runLevel on all silos, i.e. you should call it only with SKeyVal.All
    /// (at the moment there is no support for dynamic adding/removal of silos)
    void increaseRunLevel(SKey s,RunLevel newLevel){ }
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy){ }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){ }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){ }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin energy){ }
    /// a neighbor point has calculated its gradient (and energy)
    void neighborHasGradient(SKey s,LazyMPLoader!(T) p,Point[] neighbors, PointEMin energy){ }
    
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPoint(SKey s,Point,SKey owner){ }
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source,int level){ }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey,Point){ }
    // ExplorerI(T)
    /// name of this explorer
    string name(){
        return "EmptyExplorer";
    }
    /// should send an operation to evaluate to the master silos
    /// is called on all core silos in parallel
    ReturnFlag nextOp(void delegate(ExplorerI!(T)) availableAgain,int req){
        return ReturnFlag.NoOp;
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradientLocal(Point p,Real energy){
        return false;
    }
    void evaluationFailed(SKey k,Point p){ }
    
    void workOn(LocalSilosI!(T) silos){
        silos.addExplorer(this);
    }
}
