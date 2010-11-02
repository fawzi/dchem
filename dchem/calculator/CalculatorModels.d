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
import tango.io.vfs.model.Vfs;
import dchem.input.WriteOut;
import blip.parallel.mpi.MpiModels;

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
    /// a sample particle system based on real numbers (if activePrecision is Real)
    /// creation of this might be expensive, and changing it might change the internal values or not
    ParticleSys!(Real) refPSysReal();
    /// the particle system based on low precision numbers (if activePrecision is LowP)
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
    /// notification central of the current particle system & context.
    /// will receive activation, deactivation & destruction notifications
    NotificationCenter nCenter();
    /// history of the previous positions
    HistoryManager!(LowP) posHistory();
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
