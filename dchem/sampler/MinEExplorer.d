module dchem.sampler.MinEExplorer;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;

alias ulong PKey; /// key of a main point (64 bit as trade off between collision probability and being small)
alias ulong EKey; /// key of a MinEExplorer instance

struct Neighbor{
    PKey key;
    uint  direction; // uint.max means unknown
}

class MainPoint(T){
    PKey key;
    Real energies[]; /// energies of the other points (might be unallocated)
    ParticleSys!(T) pos; // position+deriv+energy
    /// 0- core, i, direction i: particle,x,y,z...
    /// flags values: 0-not explored, 1- explored by others, 2- explored close, 3- explored far
    BitArray flags;
    T[] stepScales;
    T[] repulsionScales;
    size_t nNeigh;
    Neighbor[] _neighbors;
    
    
}

// ops newpointUsed -> update neighs & flags of others,+

struct Point(T){
    MainPoint!(T) context;
    uint direction;
}

// object that keeps the journal of the computations done
class MinEJournal{
    char[] journalName;
    BinSink jSink;
    SBinSerializer jSerial;
    const char[] jVersion="MinEJournal v1.0";
    struct MainPointAdd(T){
        ulong nr;
        PKey key;
        ParticleSys!(T) pos;
        T stepSize;
        T repulsionSize;
        
        void reapply(MinEExplorer expl){}
    }
    struct EFEval(T){
        PKey key;
        ParticleSys!(T) pos;
        Real energy;

        void reapply(MinEExplorer expl){}
    }
    struct EEval(T){
        Point!(T) p;
        Real energy;

        void reapply(MinEExplorer expl){}
    }
    struct NeighAdd(T){
        Point!(T) p1;
        Point!(T) p2;

        void reapply(MinEExplorer expl){}
    }
    
}

class MinEExplorer(T): Sampler{
    char[] nameId;
    EKey key;
    MainPoint!(T)[PKey] localPoints;
    EKey[PKey] owner;
    CalculationInstance[Point] calcInProgress;
    T d;
    char[] trajDir;
    InputField evaluator;
    char[] cClass;
    
    mixin(serializeSome("dchem.MinEExplorer_"~T.stringof,
        `trajDir: directory where to store the trajectory (journal)`));
    mixin printOut!();
    void run(){
        bool restarted=false;
        evaluator.method.setupCalculatorClass();
        cClass=evaluator.method.calculatorClass;
        // possiby restarts
        if (! restarted){
            cInstance=getInstanceForClass(InstanceGetFlags.ReuseCache|InstanceGetFlags.NoAllocSubOpt|InstanceGetFlags.Wait);
        }
    }
    
    void stop(){
        
    }
    
    bool verify(CharSink s){ return true; }
}
