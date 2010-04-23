module dchem.sampler.MinEExplorer;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import Atomic=blip.sync.Atomic;
import blip.t.math.Math:max;

/// unique key for a point, dimensions might change
/// key might have some special values,
///   0: means no point,
///   1: pseudo key for points whose exploration has been excluded (for example due to constraints)
/// the key itself has the first 10 bits (size might change) that are used for collision reduction
/// and the rest that encodes a mostly sequential number (this is indirectly useful as the exploration is not random)
/// all keys with 0 at the beginning and only the first 10 bit set are thus never generated,
/// and might be used as special keys (see the special values above)
struct PointKey{
    ulong data;
    ulong key(){
        return data;
    }
    /// sets base key (should not be too large)
    void key(ulong k){
        assert(k<=0x0000_0FFF_FFFF_FFFFUL,"key too large");
        data=k&0x0000_0FFF_FFFF_FFFFUL; // drop the extra &?
    }
    /// set base key without loosing too many bits (and guranteeing that if it was non 0 it stays non 0)
    static PointKey opCall(ulong k){
        PointKey res;
        res.data=(k & 0x0000_0FFF_FFFF_FFFFUL);
        auto r=((k>>21)&0x0000_07FF_FF80_0000UL);
        if (r!=0){
            res.data=(res.data^r)|0x0000_0800_0000_0000UL;
        }
        return res;
    }
    bool isValid(){
        return (data&0x0000_0FFF_FFFF_F000UL)!=0;
    }
    mixin(serializeSome("dchem.PointKey","data"));
    mixin printOut!();
}

/// unique key for a point, and a direction in the same structure
struct PointKeyAndDir{
    /// 16 bit direction, 36 bit sequence number, 12 bit collision reduction, the bit grouping might change, don't count on it
    ulong data; /// key of a main point (64 bit as trade off between collision probability and being small). 0 is an invalid key
    static PointKeyAndDir opCall(PointKey k,uint dir){
        PointKeyAndDir res;
        if (dir!=uint.max){
            assert(dir<0xF_FFFF,"direction too large");
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)dir)<<44);
        } else {
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
        return res;
    }
    /// direction
    uint dir(){
        uint res=cast(uint)(data>>44);
        return ((res!=0xF_FFFF)?res:uint.max);
    }
    //sets direction
    void dir(uint d){
        if (d!=uint.max){
            assert(d<0xF_FFFF,"direction too large");
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)d)<<44);
        } else {
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
    }
    /// base key
    ulong key(){
        return 0x0000_0FFF_FFFF_FFFFUL & data;
    }
    /// sets base key, dropping extra bits
    void key(ulong k){
        data=(0xFFFF_F000_0000_0000 & data)|(k & 0x0000_0FFF_FFFF_FFFFUL);
    }
    /// base key
    PointKey pointKey(){
        return PointKey(0x0000_0FFF_FFFF_FFFFUL & data);
    }
    /// sets base key, dropping extra bits
    void pointKey(PointKey k){
        data=(0xFFFF_F000_0000_0000 & data)|(k.data & 0x0000_0FFF_FFFF_FFFFUL);
    }
    mixin(serializeSome("dchem.PointKeyAndDir","data"));// split in pointKey and dir???
    mixin printOut!();
}

alias ulong EKey; /// key of a MinEExplorer instance, 0 is an invalid key

// exploration:
// find smallest energy still "free"
// evaluate next direction of it
// apply constraints, if movement is too small declare it as fully visited
// bcast point as explored
// start calculation
// if too close to existing points stop calculation???
// when calculation is finished
// if mainpoint:
//    gradient -> orient neighbors, compile visited flags, perform topology analysis (attractor, min,max,...)
//    second deriv check
// else store energy, first deriv check??

/// flags about a direction
enum DirFlags{
    Free=0,
    ExploredByOthers=1,
    ExploredLocal=2,
    ExploredAndSpawn=3
}

/// an array with flags for each direction, does not use BitArray to have a fast earch of free directions
/// in large almost full arrays
class FlagsArray{
    enum{ 
        usedBits=2,
        bitsPerEl=size_t.sizeof*8/usedBits,
        mask=((~(cast(size_t)0))>>(size_t.sizeof*8-usedBits))
    }
    static if (size_t.sizeof==4){
        enum :size_t{
            add0        =0b10101010_10101010_10101010_10101011,
            bitToCheck  =0b10101010_10101010_10101010_10101010,
        }
    } else {
        enum :size_t{
            add0        =0b10101010_10101010_10101010_10101010_10101010_10101010_10101010_10101011,
            bitToCheck  =0b10101010_10101010_10101010_10101010_10101010_10101010_10101010_10101010,
        }
    }
    size_t nFlags;
    size_t[] data;
    /// construct an array with the given size
    this(size_t dim){
        nFlags=dim;
        data=new size_t[]((dim+bitsPerEl-1)/bitsPerEl);
    }
    /// generates a random FlagsArray array
    static FlagsArray randomGenerate(Rand r,ref bool acceptable){
        auto dim=generateSize(r,size_t.sizeof);
        auto res=new FlagsArray(dim);
        mkRandomArray(r,res.data,acceptable);
        return res;
    }
    // returns the flags at the given index
    uint opIndex(size_t idx){
        assert(idx<nFlags,"flags index out of bounds");
        auto block=data[idx/bitsPerEl];
        block>>=usedBits*(idx%bitsPerEl);
        return (mask&block);
    }
    /// atomic cas on a flag 
    uint atomicCAS(size_t idx,uint newF,uint oldF){
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        auto oldVal=data[blockIdx];
        assert((newF&mask)==newF,"new value is too large");
        auto localVal=(cast(size_t)newF)<<inBlockIdx;
        while(1){ // only for collisions with other updates
            auto oldFlags=(oldVal>>inBlockIdx)&0b11;
            if (oldFlags!=oldF || newF==oldF) return oldFlags;
            auto actualVal=Atomic.atomicCAS(data[blockIdx],((oldVal&inBlockMask)|localVal),oldVal);
            if (actualVal==oldVal) return oldFlags;
            oldVal=actualVal;
        }
    }
    /// atomic operation on a flag 
    void atomicOp(size_t idx,uint delegate(uint) op){
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        Atomic.atomicOp(data[blockIdx],delegate size_t(size_t oldVal){
            auto oldFlags=(oldVal>>inBlockIdx)&0b11;
            auto val=op(oldVal);
            assert((val&mask)==val,"val is too large");
            auto localVal=(cast(size_t)val)<<inBlockIdx;
            return ((oldVal&inBlockMask)|localVal);
        });
    }
    /// sets the flags at the given index (atomic, updates to *other* indexes cannot clash)
    void opIndexAssign(uint val,size_t idx){
        assert((val&mask)==val,"val is too large");
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        auto localVal=((cast(size_t)val)<<inBlockIdx);
        atomicOp(data[blockIdx],delegate size_t(size_t oldVal){
            return ((oldVal&inBlockMask)|localVal);
        });
    }
    /// returns a value different from 0 if the value val contains 00 (aligned)
    /// 4 ops + compare to 0 -> checks that 16 or 32 positions are not 00
    /// some nice bit trickery by me (fawzi :), actually the positions of the 1 bit corresponds to
    /// the start of the 00 in all cases but xxx0100 that is mapped to yyy1010 
    /// (but I did not try to take advantage of that)
    static size_t contains00(size_t val){
        auto r=((((val+add0)|val)^val)&bitToCheck);
        return r;
    }
    /// looks for a free (00) index starting at startPos.
    /// wraps around, looks at nFlags-1 as last direction, returns nFlags if no direction was found
    size_t findFreeAndSet(size_t startPos,uint newVal=0){
        if (nFlags==0) return 0;
        assert(startPos<nFlags,"startPos out of bounds");
        auto idx=startPos;
        auto blockIdx=idx/bitsPerEl;
        auto rest=idx%bitsPerEl;
        auto val=(data[blockIdx])>>(rest*usedBits);
        blockIdx+=1;
        auto last=max(nFlags-2,startPos);
        auto stopBlock=last/bitsPerEl;
        while (1){
            for (int i=rest;rest!=0;--rest){
                if ((val&0b11)==0){
                    if (atomicCAS(idx,newVal,0)==0){
                        return idx;
                    }
                }
                idx+=1;
                if (idx==startPos) { // finished search
                    if (opIndex(nFlags-1)==0){
                        return nFlags-1;
                    }
                    return nFlags;
                }
                if (idx>last) {
                    idx=0;
                    blockIdx=0;
                    stopBlock=startPos/bitsPerEl;
                    break;
                }
            }
            while(blockIdx!=stopBlock){
                val=data[blockIdx];
                auto r=((((val+add0)|val)^val)&bitToCheck);
                if (r!=0) break;
                blockIdx+=1;
            }
            idx=blockIdx*bitsPerEl;
            rest=bitsPerEl;
            val=(data[blockIdx]);
            blockIdx+=1;
        }
    }
    /// number of flags stored
    size_t length(){
        return nFlags;
    }
}

/// scale factors, put like this, so that they can be catched by the serializer
class ScaleFactors(T){
    char[] idName;
    BulkArray!(T) data;
    this(char[]idName,size_t dim){
        this.idName=idName;
        this.data=BulkArray!(T)(dim);
    }
    mixin(serializeSome("dchem.MinE.ScaleFactors("~T.stringof~")","data"));
    mixin printOut!();
}

// ops newpointUsed -> update neighs & flags of others,+
class ExpandedPoint(T){
    
}

/// a main evaluation point
class MainPoint(T){
    enum GFlags{
        EnergyOnly=(1<<0), /// only energy, no gradient (tipically results of explorations with other methods)
        Evaluated=(1<<1),  /// the energy and gradient were evaluated, frame of reference is established
        DoNotExplore=(1<<2),   /// no further exploration should be performed starting from this point
        FullyExplored=(1<<2),  /// all directions around this point have been explored
        FullyEvaluated=(1<<3), /// all directions around this point have been explored and evaluated
        OldApprox=(1<<4),      /// energy/gradient are based on old data, newer better data should be available
        LocalCopy=(1<<5),      /// this MainPoint is a local copy done for efficency reasons
    }
    MinEExplorer!(T) localContext;
    PKey key;
    Real energies[]; /// energies of the other points (might be unallocated)
    ParticleSys!(T) pos; // position+deriv+energy
    /// 0- core, i, direction i: particle,x,y,z...
    /// flags values: 0-not explored, 1- explored by others, 2- explored local, 3- explored &spawn
    FlagsArray dirFlags;
    uint gFlags;
    ScaleFactors!(T) stepScales;
    ScaleFactors!(T) repulsionScales;
    size_t nNeigh;
    Neighbor[] _neighbors;
    
    enum ExplorationProgress{
        NoDirAvailable,
        WaitingEvaluationCompletion,
        Explorating,
    }
    /// explores the next direction, immediately marks it as in exploration
    PointKeyAndDir exploreNext(bool cheapGrad){
        if ((gFlags&GFlags.DoNotExplore)!=0){
            return SubPoint;
        }
        if ((gFlags&GFlags.Evaluated)==0){
            if (!localContext.isInProgress(Point(PKey,0))){
                synchronized(this){
                    if ((gFlags&GFlags.Evaluated)==0){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("non evaluated ans non in progress ")(key)(" in ")(localContext.nameId);
                        }));
                    }
                }
            }
            return false; // needs to wait evaluation
        }
        if ((gFlags&GFlags.FullyExplored)==0){
            synchronized(this){
                if ((gFlags&GFlags.FullyExplored)==0){
                    if(dirFlags.atomicCAS(1,DirFlags.None,DirFlags.ExploredLocal)==DirFlags.None){
                        //exploreDown
                    }
                }
            }
        }
        return PointKeyAndDir(0,0);
    }
    /// spawn exploration in the given direction
    MainPoint spawnDirection(uint dir){
        return false;
    }
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
        PSysWriter!(T) pos;
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
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    EKey key;
    MainPoint!(T)[PKey] localPoints;
    EKey[PKey] owner;
    CalculationInstance[Point] calcInProgress;
    Random randLocal;
    RandomSync rand;
    
    T d;
    char[] trajDir;
    InputField evaluator;
    char[] cClass;
    /// returns a globally unique string 
    char[] nextUniqueStr(){
        return collectAppender(delegate void(CharSinker s){
            s(nameId); s("_"); writeOut(s,nextUniqueId.next());
        });
    }
    /// returns a most likely valid point id
    PKey nextPointId(){
        ushort r;
        rand(r);
        return ((nextPntNr.next())<<16)|cast(ulong)r;
    }
    /// returns true if the evaluation of the given point is in progress
    bool isInProgress(Point!(T)p){
        bool res=false;
        synchronized(this){
            res=(p in calcInPyrogressy)!is null;
        }
        return this;
    }
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
