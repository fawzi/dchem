module dchem.calculator.CalculatorModels;
import dchem.input.RootInput;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import blip.serialization.Serialization;
import blip.text.TextParser;
import blip.io.BasicIO;
import blip.io.Console;
import dchem.input.ReadIn;
import dchem.Common;
import blip.util.TangoLog;
import blip.core.Variant;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import tango.io.FilePath;
import dchem.input.ReadIn2PSys;
import blip.util.NotificationCenter;
import dchem.sys.PIndexes;
import blip.BasicModels;
import dchem.sys.Constraints;
import dchem.sys.DynVars;
import tango.io.vfs.model.Vfs;
import dchem.input.WriteOut;
import blip.parallel.mpi.MpiModels;
import blip.container.BulkArray;

/// a task that can be excuted locally on the calculation context
interface RemoteCCTask: Serializable{
    /// starts the task with the given LocalCalculationContext
    void workOn(LocalCalculationContext);
    /// might stop the task, or not, might return immediatly (even if the task is still running)
    void stop();
}

/// calculator setup (chooses the method to perform each calculation)
interface Method:InputElement{
    /// activates the method in the context of the given parallel environment
    /// this can be used to do setups shared by all contexts of this method
    /// it might be called several times, but should always have the same arguments
    void setup(LinearComm pEnv, CharSink log);
    /// gets a calculator to perform calculations with this method, if possible reusing the given history
    /// if wait is true waits until a context is available
    CalculationContext getCalculator(bool wait,ubyte[]history);
    /// drops the history with the given id
    void dropHistory(ubyte[]history);
    /// clears all history
    void clearHistory();
}

/// amount of change since the last calculation in the context
enum ChangeLevel{
    FirstTime=0,
    AllChanged=1,
    PosChanged=2,
    SmallPosChange=3,
    SmoothPosChange=4
}

enum Precision{
    Real,
    LowP
}

/// input object that generates the objects that defines the macroscopic distance space (DistOps)
interface DistOpsGen:InputElement{
    DistOps DistOpsForContext(CalculationContext);
}
/// an object that handles the various distance measures, and distance related things
/// it defines the macroscopic distance space, the microscopic distance space (for small changes)
/// is defined in ParticleSys by its overlap, here one defines things like PBC
interface DistOps{
    /// wraps deltaX (a difference between two X points) so that it is as small as possible
    /// compatibly with the implicit symmetries (but not with the explicit ones)
    void wrap(ParticleSys!(Real)pSys,DynPVector!(Real,XType)deltaX);
    /// ditto
    void wrap(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)deltaX);
    /// ditto
    void wrapReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)deltaX);
    /// ditto
    void wrapLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)deltaX);
    
    /// get the periodic copy of x that is closest (in first image convetion, i.e. in the h_inv distance)
    /// to pSys (this might be different than the closest one for very skewed h)
    void makeClose(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x);
    /// ditto
    void makeClose(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x);
    /// ditto
    void makeCloseReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x);
    /// ditto
    void makeCloseLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x);
    
    /// distance in the reduced units, this is the norm2 distance between x and pSys.dynVars.x after a makeClose
    /// or deltaX.norm2 after calling wrap, where deltaX=x-pSys.dynVars.x
    /// distances larger than threshold might not be computed accurately
    Real reducedDist(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold=Real.max);
    /// ditto
    Real reducedDist(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold=Real.max);
    /// ditto
    Real reducedDistReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold=Real.max);
    /// ditto
    Real reducedDistLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold=Real.max);
    
    /// distance in reduced units between *one* particle of kinf k at pos1,ord1,dof1, and possibly many
    /// particles of the same kind at pos2,ord2,dof2
    void rDistOneToN(ParticleSys!(Real)pSys,KindIdx k,
        BulkArray!(Vector!(Real,3))pos1,BulkArray!(Quaternion!(Real))ord1,BulkArray!(Real)dof1,
        BulkArray!(Vector!(Real,3))pos2,BulkArray!(Quaternion!(Real))ord2,BulkArray!(Real)dof2,
        Real[] dists);
    /// ditto
    void rDistOneToN(ParticleSys!(LowP)pSys,KindIdx k,
        BulkArray!(Vector!(LowP,3))pos1,BulkArray!(Quaternion!(LowP))ord1,BulkArray!(LowP)dof1,
        BulkArray!(Vector!(LowP,3))pos2,BulkArray!(Quaternion!(LowP))ord2,BulkArray!(LowP)dof2,
        LowP[] dists);
    /// ditto
    void rDistOneToNReal(ParticleSys!(Real)pSys,KindIdx k,
        BulkArray!(Vector!(Real,3))pos1,BulkArray!(Quaternion!(Real))ord1,BulkArray!(Real)dof1,
        BulkArray!(Vector!(Real,3))pos2,BulkArray!(Quaternion!(Real))ord2,BulkArray!(Real)dof2,
        Real[] dists);
    /// ditto
    void rDistOneToNLowP(ParticleSys!(LowP)pSys,KindIdx k,
        BulkArray!(Vector!(LowP,3))pos1,BulkArray!(Quaternion!(LowP))ord1,BulkArray!(LowP)dof1,
        BulkArray!(Vector!(LowP,3))pos2,BulkArray!(Quaternion!(LowP))ord2,BulkArray!(LowP)dof2,
        LowP[] dists);
    
    /// full (cartesian) distance, this is the reducedDist in the full system using cartesian units
    /// in general it might be expensive to calculate, or be not better than the 
    /// if threshold is different from 0, then as soon as the distance is detected to be larger than it
    /// it is returned, which might be quicker if one does not care when for the exact value of distances
    /// larger than threshold
    Real fullDist(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold=Real.max);
    /// ditto
    Real fullDist(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold=Real.max);
    /// ditto
    Real fullDistReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x, Real threshold=Real.max);
    /// ditto
    Real fullDistLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x, Real threshold=Real.max);
}

/// input object that generates the explicit symmetries (SymmNeighLoopers)
interface SymmNeighLooperGen:InputElement{
    SymmNeighLooper symmNeighLooperForContext(CalculationContext);
}
/// an object that loops on all explicit symmetry equivalent structures generated by neigh that
/// are within epsilon of pSys (the exact meaning of epsilon is up to the specific operation)
interface SymmNeighLooper{
    /// loops on all explicit symmetry equivalent structures generated by neigh that are within epsilon of pSys
    int loopOnNeighWithin(ParticleSys!(Real)pSys,DistOps dOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody);
    /// ditto
    int loopOnNeighWithin(ParticleSys!(LowP)pSys,DistOps dOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody);
    /// ditto
    int loopOnNeighWithinReal(ParticleSys!(Real)pSys,DistOps dOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody);
    /// ditto
    int loopOnNeighWithinLowP(ParticleSys!(LowP)pSys,DistOps dOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody);
}
/// helper structure to support nicer foreach looping
struct SymmNeighLooperHelper(T){
    SymmNeighLooper looper;
    ParticleSys!(T)pSys;
    DynamicsVars!(T)neigh;
    T epsilon;
    int opApply(int delegate(ref DynPVector!(T,XType))loopBody){
        return looper.loopOnNeighWithin(pSys,neigh,epsilon,loopBody);
    }
}
/// function to support nicer foreach syntax for loop on symmetry equivalent neighbors
SymmNeighLooperHelper!(T) symmNeighLooper(T,U)(SymmNeighLooper looper,ParticleSys!(T)pSys,DynPVector!(T,XType)neigh,U epsilon){
    SymmNeighLooperHelper!(T) lHelper;
    lHelper.looper=looper;
    lHelper.pSys=pSys;
    lHelper.neigh=neigh;
    lHelper.epsilon=cast(T)epsilon;
    return lHelper;
}

/// calculation that might have been aready partially setup, in particular the
/// number of elements,... cannot change
/// these method can be called from a remote computer, so no overload, and thus
/// the non D like *Set methods
interface CalculationContext{
    /// the active precision, this says which one of pSysReal/pSysLowP and constraintsReal/constraintsLowP
    /// is active. The other will most likely be null.
    Precision activePrecision();
    /// unique identifier for this context
    char[] contextId();
    /// the particle system based on real numbers (if activePrecision is Real)
    PSysWriter!(Real) pSysWriterReal();
    /// the particle system based on low precision numbers (if activePrecision is LowP)
    PSysWriter!(LowP) pSysWriterLowP();
    /// the particle system based on real numbers (if activePrecision is Real)
    void pSysWriterRealSet(PSysWriter!(Real));
    /// the particle system based on low precision numbers (if activePrecision is LowP)
    void pSysWriterLowPSet(PSysWriter!(LowP));
    /// a sample particle system based on real numbers
    /// creation of this might be expensive, and it might be empty or not update the real PSys
    ParticleSys!(Real) refPSysReal();
    /// the particle system based on low precision numbers
    /// creation of this might be expensive, and it might be empty or not update the real PSys
    ParticleSys!(LowP) refPSysLowP();
    /// the generator for the constraints of this context
    ConstraintGen constraintGen();
    /// the system struct
    SysStruct sysStruct();
    /// change level since the last calculation
    ChangeLevel changeLevel();
    /// sets the change level
    void changeLevelSet(ChangeLevel);
    /// decreases the change level to at least changeLevel
    void changedDynVars(ChangeLevel changeLevel,Real diff);
    /// the total potential energy
    Real potentialEnergy();
    /// changes the position (utility method)
    void posSet(SegmentedArray!(Vector!(Real,3)) newPos);
    /// gets the positions (utility method)
    SegmentedArray!(Vector!(Real,3)) pos();
    /// changes the velocities (utility method)
    void dposSet(SegmentedArray!(Vector!(Real,3)) newDpos);
    /// gets the velocities (utility method)
    SegmentedArray!(Vector!(Real,3)) dpos();
    /// changes the forces (utility method)
    void mddposSet(SegmentedArray!(Vector!(Real,3)) newDdpos);
    /// gets the forces (utility method)
    SegmentedArray!(Vector!(Real,3)) mddpos();
    /// updates energy and/or forces
    void updateEF(bool updateE=true,bool updateF=true);
    /// called automatically after creation, but before any energy evaluation
    /// should be called before working again with a deactivated calculator
    void activate();
    /// call this to possibly compress the context and empty caches (i.e. before a long pause in the calculation)
    void deactivate();
    /// call this to gives back the context (after all calculations with this are finished) might delete or reuse it
    void giveBack();
    /// tries to stop a calculation in progress, recovery after this is not possible (only giveBack can be called)
    void stop();
    /// method that has generated this context
    Method method();
    /// stores the history somewhere and returns an id to possibly recover that history at a later point
    /// this is just an optimization, it does not have to do anything. If implemented then
    /// method.getCalculator, .dropHistory and .clearHistory have to be implemented accordingly
    ubyte[]storeHistory();
    /// url to access this from other processes
    char[] exportedUrl();
    /// execute the operation t locally on the context
    void executeLocally(RemoteCCTask t);
    /// log for the calculator context
    void logMsg1(char[]);
}

ParticleSys!(T) refPSysT(T)(CalculationContext ctx){
    static if(is(T==LowP)){
        return ctx.refPSysLowP();
    } else static if (is(T==Real)){
        return ctx.refPSysReal();
    } else {
        static assert(0,"unexpected type "~T.stringof);
    }
}

/// represent the local view to a calculation that might have been aready partially setup, in particular the
/// number of elements,... cannot change, gives direct access to modifiable ParticleSys, and more
interface LocalCalculationContext: CalculationContext{
    /// the particle system based on real numbers (if activePrecision is Real)
    ParticleSys!(Real) pSysReal();
    /// the particle system based on low precision numbers (if activePrecision is LowP)
    ParticleSys!(LowP) pSysLowP();
    /// the constraints of this context (if activePrecision is Real)
    ConstraintI!(Real) constraintsReal();
    /// the constraints of this context (if activePrecision is LowP)
    ConstraintI!(LowP) constraintsLowP();
    /// the symmNeighLooper for this context
    SymmNeighLooper symmNeighLooper();
    /// the distance operations for this context
    DistOps distOps();
    /// notification central of the current particle system & context.
    /// will receive activation, deactivation & destruction notifications
    NotificationCenter nCenter();
    /// history of the previous positions
    HistoryManager!(LowP) posHistory();
    /// the logger of this context
    CharSink logger();
}

/// templatized way to extract the pSys from a LocalCalculationContext
ParticleSys!(T) pSysT(T)(LocalCalculationContext c){
    static if (is(Real==T)){
        return c.pSysReal;
    } else static if(is(LowP==T)){
        return c.pSysLowP;
    } else {
        static assert("requested a pSys "~T.stringof~" which is neither Real ("~
            Real.stringof~") nor LowP ("~LowP.stringof~")");
    }
}
/// templatized way to extract the pSysWriter from a CalculationContext
PSysWriter!(T) pSysWriterT(T)(CalculationContext c){
    static if (is(Real==T)){
        return c.pSysWriterReal;
    } else static if(is(LowP==T)){
        return c.pSysWriterLowP;
    } else {
        static assert("requested a pSysWriterLowP "~T.stringof~" which is neither Real ("~
            Real.stringof~") nor LowP ("~LowP.stringof~")");
    }
}
/// templatized way to extract the pSys from a LocalCalculationContext
void pSysWriterSetT(T)(CalculationContext c,PSysWriter!(T) sys){
    static if (is(Real==T)){
        return c.pSysWriterRealSet(sys);
    } else static if(is(LowP==T)){
        return c.pSysWriterLowPSet(sys);
    } else {
        static assert("requested to set pSysWriter "~T.stringof~" which is neither Real ("~
            Real.stringof~") nor LowP ("~LowP.stringof~")");
    }
}
