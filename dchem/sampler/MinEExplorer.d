module dchem.sampler.MinEExplorer;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import Atomic=blip.sync.Atomic;
import blip.t.math.Math:max;
import tango.util.container.more.Heap;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;

/// unique key for a point, dimensions might change
/// key might have some special values,
///   0: means no point,
///   1: pseudo key for points whose exploration has been excluded (for example due to constraints)
/// the key itself has the first 10 bits (size might change) that are used for collision reduction
/// and the rest that encodes a mostly sequential number (this is indirectly useful as the exploration is not random),
/// and to allow efficient storage of dense arrays (making unbundling of properties of points easy)
/// all keys with 0 at the beginning and only the first 10 bit set are thus never generated,
/// and might be used as special keys (see the special values above)
struct Point{
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
    static Point opCall(ulong k){
        Point res;
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
    mixin(serializeSome("dchem.Point","data"));
    mixin printOut!();
    hash_t toHash(){
        static if(hash_t.sizeof==4){
            return cast(hash_t)(0xFFFF_FFFF&(data^(data>>32)));
        } else {
            return cast(hash_t)data;
        }
    } 
}

/// unique key for a point, and a direction in the same structure
struct PointAndDir{
    /// 16 bit direction, 36 bit sequence number, 12 bit collision reduction, the bit grouping might change, don't count on it
    ulong data; /// key of a main point (64 bit as trade off between collision probability and being small). 0 is an invalid key
    static PointAndDir opCall(Point k,uint dir){
        PointAndDir res;
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
    Point point(){
        return Point(0x0000_0FFF_FFFF_FFFFUL & data);
    }
    /// sets base key, dropping extra bits
    void point(Point k){
        data=(0xFFFF_F000_0000_0000 & data)|(k.data & 0x0000_0FFF_FFFF_FFFFUL);
    }
    mixin(serializeSome("dchem.PointAndDir","data"));// split in Point and dir???
    mixin printOut!();
    struct ExpandedPointAndDir{
        PointAndDir pDir;
        uint dir(){
            return pDir.dir();
        }
        void dir(uint d){
            pDir.dir=d;
        }
        Point point(){
            return pDir.point();
        }
        void point(Point p){
            pDir.point=p;
        }
        mixin(serializeSome("dchem.PointDir","point|dir"));// split in Point and dir???
        mixin printOut!();
    }
    ExpandedPointAndDir expanded(){
        ExpandedPointAndDir e;
        e.pDir=*this;
        return e;
    }
}

alias ulong EKey; /// key of a MinEExplorer instance, 0 is an invalid key

/// flags about a direction
enum DirFlags{
    Free=0,
    ExploredByOthers=1,
    ExploredLocal=2,
    ExploredAndSpawn=3
}

/// prints the bit pattern of a, at least minBits are printed (minBits can be at most 64)
char[] printBits(size_t a,int minBits=1){
    char[64] res;
    auto i=64;
    while(a!=0){
        --i;
        if ((a&1)==1){
            res[i]='1';
        } else {
            res[i]='0';
        }
        a>>=1;
    }
    auto mBit=((minBits<=64)?minBits:64);
    while (64-i<mBit){
        res[--i]='0';
    }
    return res[i..$].dup;
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
    this(){}
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
    uint opIndex(size_t idx)
    in{
        if (idx>=nFlags){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("flags index ")(idx)(" out of bounds (0,")(nFlags)(")");
            }),__FILE__,__LINE__);
        }
    } body {
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
            auto oldFlags=((oldVal>>inBlockIdx)&0b11);
            if (oldFlags!=oldF) return oldFlags;
            auto actualVal=Atomic.atomicCAS(data[blockIdx],((oldVal&inBlockMask)|localVal),oldVal);
            if ((actualVal&inBlockMask)==(oldVal&inBlockMask)) return ((actualVal>>inBlockIdx)&0b11);
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
        if (nFlags==1) {
            if (atomicCAS(0,newVal,0)==0){
                return 0;
            }
            return 1;
        }
        assert(startPos<nFlags,"startPos out of bounds");
        auto idx=startPos;
        auto blockIdx=idx/bitsPerEl;
        auto rest=idx%bitsPerEl;
        auto val=(data[blockIdx])>>(rest*usedBits);
        rest=bitsPerEl-rest;
        blockIdx+=1;
        auto last=((nFlags>startPos+2)?(nFlags-2):startPos);
        auto stopBlock=last/bitsPerEl;
        bool wasReset=false;
        while (1){
            while (rest!=0){
                --rest;
                if (atomicCAS(idx,newVal,0)==0){
                    return idx;
                }
                idx+=1;
                val>>=usedBits;
                if (idx==startPos) { // finished search¨
                    if (opIndex(nFlags-1)==0){
                        if (atomicCAS(nFlags-1,newVal,0)==0){
                            return nFlags-1;
                        }
                    }
                    return nFlags;
                }
                if (idx>last) {
                    if (wasReset) throw new Exception("reset loop",__FILE__,__LINE__);
                    wasReset=true;
                    idx=0;
                    blockIdx=0;
                    stopBlock=startPos/bitsPerEl;
                    if (idx==startPos) { // finished search
                        if (atomicCAS(nFlags-1,newVal,0)==0){
                            return nFlags-1;
                        }
                        return nFlags;
                    }
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
    mixin(serializeSome("dchem.FlagsArray","nFlags|data"));
    void desc(CharSink sink){
        auto s=dumper(sink);
        s("{class:\"dchem.FlagsArray\",length:")(nFlags)(",data:");
        auto i=nFlags;
        while(i!=0){
            if(i!=nFlags) s(".");
            --i;
            s(printBits(opIndex(i),2));
        }
        s("}");
    }
}

/// scale factors, put like this, so that they can be catched by the serializer, * to remove *
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

struct LocalPoint(T){
    Point point;
    PointAndDir generator;
    MainPoint!(T) _mainPoint;
    MainPoint!(T) mainPoint(){
        if (_mainPoint !is null) return _mainPoint;
        synchronized{
            if (_mainPoint !is null) return _mainPoint;
            if (! generator.isValid()) throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)(point)(" has no mainPoint and generator ")(generator)(" is non valid");
            }),__FILE__,__LINE__);
            
        }
    }
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
    Point point;
    /// position of the point (energy, derivatives,... might be invalid, only real positions have to be valid)
    ParticleSys!(T) pos;
    /// direction of the minimum in the dual space with norm 1 wrt. euclidische norm
    /// (frame of reference for minimization)
    DynPVect!(T,2) minDir;
    /// directions, at index 0- core, index i, direction i of the derivatives in the dual space (DynPVect!(T,2))
    /// flags values are the DirFlags: 0-not explored, 1- explored by others, 2- explored local, 3- explored &spawn
    FlagsArray dirFlags;
    /// localNeighbors, point is the evaluation point, dir the direction with respect to this point to generate it
    /// energy its energy
    GrowableArray!(PointDirAndEnergy) localNeighbors;
    /// external neighbors, and the direction with respect to this point in which they are
    GrowableArray!(PointAndDir) externalNeighbors;
    /// energy of this point
    Real energy0;
    /// scale of mindir to recover the dual gradient
    T minDirScale;
    /// repulsion size (useful to block out regions)
    T repulsionSize;
    /// exploration size for this point (useful to recover positions) 
    T explorationSize;
    /// bit or of GFlags of the current point
    uint gFlags;
    
    /// explores the next direction, immediately marks it as in exploration
    /// returns an invalid point only if all direction are explored, and an invalid direction but this point if a
    /// frame of exploration is not yet defined
    PointAndDir exploreNext(bool cheapGrad){
        auto exploreFlags=DirFlags.ExploredLocal;
        if (cheapGrad) exploreFlags=DirFlags.ExploredAndSpawn;
        if ((gFlags&GFlags.DoNotExplore)!=0){
            return PointAndDir(0,uint.max);
        }
        if ((gFlags&GFlags.Evaluated)==0){
            if (!localContext.isInProgress(PointAndDir(point,0))){
                synchronized(this){
                    if ((gFlags&GFlags.Evaluated)==0){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("non evaluated and non in progress ")(point)(" in ")(localContext.nameId);
                        }));
                    }
                }
            }
            return PointAndDir(this.point,uint.max); // needs to wait evaluation
        }
        if ((gFlags&GFlags.FullyExplored)==0){
            synchronized(this){
                if ((gFlags&GFlags.FullyExplored)==0){
                    size_t nextDir=dirFlags.length;
                    if(dirFlags.findFreeAndSet(1,exploreFlags)){
                        nextDir=0;
                    } else {
                        auto nextDir=dirFlags.findFreeAndSet(start,exploreFlags);
                    }
                    if (nextDir!=dirFlags.length) return PointAndDir(this.point,nextDir);
                }
            }
        }
        return PointAndDir(0,0);
    }
    
    /// calculates the point in the given direction
    ParticleSys!(T) posInDir(uint dir){
        if (dir==0) return pos;
        rDir=(dir-1)/2;
        T val=1;
        if ((dir&1)==0) val=-1;
        if (dir==dirFlags.length-1){
            rDir=0;
            val=-1;
        }
        auto v0=minDir.emptyCopy;
        v0[]=0;
        v0[dir]=1;
        auto newDir=rotateEiV(V,M)(0,minDir,v0);
        newDir*=explorationSize;
        auto resPSys=pos.dup(PSDupLevel.DynProperties|PSDupLevel.DynPNullDx|PSDupLevel.DynPNullMddx|PSDupLevel.HiddenVars);
        void addFromDualTSpace(T)(newDir,resPSys.dynVars.x);
        resPSys.updateHVars();
        return null;
    }
}

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
        Point!(T) p2;

        void reapply(MinEExplorer expl){}
    }
    
}

class MinHeapSync(T){
    MinHeap!(T) heap;
    WaitCondition nonEmpty;
    
    bool nonEmptyHeap(){
        return heap.length!=0;
    }
    
    this(){
        nonEmpty=new WaitCondition(&nonEmptyHeap);
    }
    void push(T[] t){
        synchronized(this){
            heap.push(t);
        }
        nonEmpty.checkCondition();
    }
    void push(T t){
        synchronized(this){
            heap.push(t);
        }
        nonEmpty.checkCondition();
    }
    T pop(){
        synchronized(this){
            return heap.pop();
        }
    }
    /// returns the minimal energy elements, waits if no elements is available until some becomes available
    T waitPop(){
        while (1){
            synchronized(this){
                if (heap.length>0)
                    return heap.pop();
            }
            nonEmpty.wait();
        }
    }
}

struct PointDirAndEnergy{
    Real energy;
    PointAndDir point;
    mixin(serializeSome("dchem.PointDirAndEnergy","energy|point"));
}

interface ActiveExplorer(T){
    /// stops all explorers (at the moment there is no support for dynamic adding/removal of explorers, worth adding???)
    void shutdown(EKey);

    /// key of this explorer
    EKey key();
    /// where the journal is kept
    char[] journalPos();
    
    /// adds energy to the evaluated ones (and removes from the inProgress)
    void addEnergyEval(EKey,PointAndEnergy p);
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void addExploredPoint(EKey,Point key,PSysWriter!(T) pos,uint flags);
    /// finished exploring the given point
    void finishedExploringPoint(EKey,Point);
    /// drops all calculation/storage connected with the given point (called upon collisions)
    void dropPoint(EKey,Point);
    
    /// returns the position of the given point
    /// Pos
}
alias HashSet Set;
enum { batchSize=512 }
/// dense array of elements indexed by local point keys
class DenseLocalPointArray(T){
    BatchedGrowableArray!(Point,batchSize) keys;
    BatchedGrowableArray!(T,batchSize) values;
    
    T opIndex(Point k){
        synchronized(this){
            auto idx=lbound(keys,k);
            assert(idx<keys.length && keys[idx]==k,"could not find local key "~k.toString());
            if (values.length<idx){
                return T.init;
            }
            return values[idx];
        }
    }
    T* ptrI(Point k){
        auto idx=lbound(keys,k);
        assert(idx<keys.length && keys[idx]==k,"could not find local key "~k.toString());
        if (values.length<idx){
            values.growTo(keys.length);
        }
        return &(values[idx]);
    }
    void opIndexAssign(T val,Point k){
        synchronized(this){
            *ptrI()=val;
        }
    }
    /+static if (T U:U*){
        int opApply(int delegate(ref U el)){
            synchronized(this){
                auto k=keys.view();
                values.growTo(k.length);
                auto v=values.view();
            }
            //v.sLoop()
        }
    }+/
}

class MinEExplorer(T): Sampler,ActiveExplorer!(T){
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    EKey key;
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    /// list of all the currently active explorers (including this one)
    ActiveExplorer!(T)[] activeExplorers;
    BatchGrowableArray!(Point,batchSize) localPointsKeys; // indexes for local owned points
    HashMap!(Point,EKey) owner;
    CalculationInstance[Point] localCalcInProgress;
    Set!(Point) calcInProgress;
    Random randLocal;
    RandomSync rand;
    /// step used for the discretization
    T discretizationStep;
    /// place to store the trajectory (journal)
    char[] trajDir;
    /// the current method for energy evaluation
    InputField evaluator;
    
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
        bool restarted=false;
        evaluator.method.setupCalculatorClass();
        // possiby restarts
        if (! restarted){
            cInstance=getInstanceForClass(InstanceGetFlags.ReuseCache|InstanceGetFlags.NoAllocSubOpt|InstanceGetFlags.Wait);
        }
        master.nextCalculation();
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
        
    }
    void stop(){
        
    }
    
    bool verify(CharSink s){ return true; }
}
