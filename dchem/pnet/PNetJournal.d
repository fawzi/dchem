/// keeps a journal of the changes
module dchem.pnet.PNetJournal;

// object that keeps the journal of the computations done
class MinEJournal{
    char[] journalName;
    BinSink jSink;
    SBinSerializer jSerial;
    const char[] jVersion="MinEJournal v1.0";
    struct MainPointAdd(T){
        Point point;
        PSysWriter!(T) pos;
        T stepSize;
        T repulsionSize;
        
        void reapply(MinEExplorer expl){}
    }
    struct EFEval(T){
        Point key;
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
        Point!(T)[] p2;

        void reapply(MinEExplorer expl){}
    }
    
}
