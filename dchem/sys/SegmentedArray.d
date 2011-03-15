/// segmented array: efficient dense storage for particle properties that have kind dependent storage need
/// this is the basis of for storing most properties, and is used for position,charge,... arrays
///
/// as extra optimization one could allocate the various struct arrays at the same time (giving a "nicer" memory layout)
module dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import blip.container.BulkArray;
import dchem.sys.SubMapping;
import dchem.Common;
import blip.narray.NArray;
import blip.sync.Atomic;
import blip.core.Traits;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import blip.core.Variant;
import blip.parallel.smp.WorkManager;
import blip.container.AtomicSLink;
import blip.container.Pool;
import blip.io.Console;
import blip.util.Convert;
import blip.core.Traits;
import blip.container.Cache;
import blip.stdc.string: memcpy;
import blip.math.Math: max;

enum ParaFlags{
    FullPara,
    KindPara,
    Sequential,
    DataLoop=1<<5,
    MaskBase= ~DataLoop,
}

/// structure of a segmented array, kindStarts is valid only after freezing
final class SegmentedArrayStruct{
    enum Flags{
        None=0,   /// no special flags
        Frozen=1, /// kind dimensions cannot be changed anymore (this is set when an array is created)
        Direct=2, /// =submapping.mappingKind & MappingKind.Direct, worth caching?
        Min1=4,   /// at least one value per kind should be stored
    }
    char[] name; /// a name (for debugging)
    SubMapping submapping;
    KindRange  kRange;
    index_type[]   _kindDims; /// number of basic elements for each particle of this kind
    index_type[]   kindStarts; /// starts of the kinds
    index_type[]   kindIncrements; /// when looping the increment to add for each particle of the given kind
    index_type[]   mKindDims; /// minimum number of basic elements (this is max(1,_kindDims) if Min1, otherwise _kindDims)
    Flags flags;
    mixin(serializeSome("dchem.sys.SegmentedArrayStruct","Describes the structure of a SegmentedArray",
        "submapping|kRange|_kindDims|kindStarts"));
    mixin printOut!();
    /// allocates a new SegmentedArray with the given kind dimensions
    /// if the Min1 flag is set (default) then at least one value per kind is stored 
    this(char[]name,SubMapping submapping,KindRange kRange,index_type[]kindDims=null,Flags f=Flags.Min1){
        assert(submapping !is null,"submapping needed to allocate");
        this.name=name;
        this.submapping = submapping;
        this.kRange     = kRange    ;
        assert(kRange!=KindRange.all,"the range should be explicitly set (to avoid missing empty kinds)");
        assert(submapping.mappingKind&MappingKind.KindPreserving,"submapping must be KindPreserving");
        this.flags=(f& ~(Flags.Direct|Flags.Frozen))|(((submapping.mappingKind & MappingKind.Direct)!=0)?
                Flags.Direct:Flags.None);
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        if (kindDims.length==0 && nkinds!=0){
            kindDims=new index_type[](nkinds);
            kindDims[]=0; // should be the default
        }
        this._kindDims   = kindDims  ;
        this.kindIncrements=new index_type[](nkinds);
        this.mKindDims=new index_type[](nkinds);
        assert(_kindDims.length==nkinds,"kindDims has the wrong length");
        assert(kRange in submapping.lKRange,"submapping is smaller than current kRange");
        //if ((flags & Flags.Frozen)!=0){
        //    recalculateStarts();
        //} else {
            kindStarts[]=index_type.max;
        //}
    }
    /// number of particles for the given kind
    /// make this faster???
    index_type nParticles(KindIdx k){
        assert(k in kRange);
        auto kLoc=k-submapping.lKRange.kStart;
        return submapping.kindStarts[kLoc+1]-submapping.kindStarts[kLoc];
    }
    /// ditto
    this(char[] name,SubMapping submapping, KindRange kRange, NArray!(index_type,1) kindDims,Flags f=Flags.Min1){
        if (kindDims is null) {
            this(name,submapping,kRange,cast(index_type[])null,f);
        } else {
            index_type[] kDims=new index_type[kindDims.shape[0]];
            for (size_t i=0;i<kindDims.shape[0];++i){
                kDims[i]=kindDims[i];
            }
            this(name,submapping,kRange,kDims,f);
        }
    }
    // internal for serialization
    this(){ }
    /// makes the structure non modifiable (done before creating data arrays)
    typeof(this) freeze(){
        if ((flags & Flags.Frozen)==0){
            synchronized(this){
                if ((flags & Flags.Frozen)==0){
                    recalculateStarts();
                    flags|=Flags.Frozen;
                }
            }
        }
        return this;
    }
    /// recalculates the starts of the array
    void recalculateStarts(){
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        kindStarts.length=nkinds+1;
        
        kindStarts[0]=0;
        bool min1=((flags & Flags.Min1)!=0);
        auto kindShift=cast(size_t)(kRange.kStart-submapping.lKRange.kStart);
        for (size_t i=0;i<nkinds;++i){
            if (min1 && _kindDims[i]==0){
                kindStarts[i+1]=kindStarts[i]+1; // always keep a "kind owned" value
                kindIncrements[i]=0;
                mKindDims[i]=1;
            } else {
                kindStarts[i+1]=kindStarts[i]
                    +_kindDims[i]*(submapping.kindStarts[i+1+kindShift]-submapping.kindStarts[i+kindShift]);
                kindIncrements[i]=_kindDims[i];
                mKindDims[i]=_kindDims[i];
            }
        }
    }
    /// equality, add also equivalence and inclusion??
    equals_t opEquals(Object o2){
        SegmentedArrayStruct s2=cast(SegmentedArrayStruct)o2;
        if (s2!is null){
            return this is s2 || (kRange==s2.kRange && kindStarts==s2.kindStarts
                && submapping==s2.submapping && flags==s2.flags);
        }
        return false;
    }
    /// return the dimension of kind i
    index_type kindDim(T)(T k){
        assert(cast(KindIdx)k in kRange,"kind out of range");
        return _kindDims[cast(size_t)(k-kRange.kStart)];
    }
    /// sets the dimension of kind i, dangerous when called concurrently
    void setKindDim(KindIdx k,index_type val){
        if((flags & Flags.Frozen)==0)
            throw new Exception("index assign on frozen struct",__FILE__,__LINE__);
        assert(k in kRange,"kind out of range");
        auto i=cast(size_t)(k-kRange.kStart);
        return _kindDims[i]=val;
    }
    /// adds the the dimension, returns the original value
    index_type addToKindDim(KindIdx k,index_type val){
        assert(k in kRange,"kind out of range");
        auto i=cast(size_t)(k-kRange.kStart);
        if((flags & Flags.Frozen)!=0)
            throw new Exception("index assign on frozen struct",__FILE__,__LINE__);
        return atomicAdd(_kindDims[cast(size_t)i],val);
    }
    /// returns a non frozen copy
    typeof(this) dup(){
        return new SegmentedArrayStruct(name~"Dup",submapping,kRange,_kindDims.dup,flags & (~Flags.Frozen));
    }
    /// the length of the data stored in the segmented array (i.e. sArr.data.length), change???
    size_t dataLength(){
        return kindStarts[$-1]-kindStarts[0];
    }
    /// length of a particle based loop on this (i.e skipping kinds without particles, and counting
    /// the number of particles for dimension of size 0 if Min1, i.e. a single value per kind)
    /// this is the size of sLoop and pLoop. The real storage size is .data.length
    size_t length(){
        size_t iterRef=this.dataLength;
        auto submap=submapping;
        if((this.flags&SegmentedArrayStruct.Flags.Min1)!=0){
            for (auto k=this.kRange.kStart;k<this.kRange.kEnd;++k){
                if (this.kindDim(k)==0){
                    iterRef-=1;
                    iterRef+=submap.kindStarts[k-submap.lKRange.kStart+1]-submap.kindStarts[k-submap.lKRange.kStart];
                }
            }
        }
        return iterRef;
    }
}

/// maps each particle kind to a contiguous memory segment
/// useful to put several segmentedArrays on the same memory, sharing full system and pieces of it
class SegArrMemMap(T):PoolI!(SegmentedArray!(T)){
    SegmentedArrayStruct arrayStruct;
    size_t[] kindOffsets; /// start of each kind wrt. to the guard start (in bytes)
    size_t[] kindByteIncrements; /// increments in bytes
    size_t allocSize;     /// alloc size of the BulkArray (in bytes)
    size_t alignBytes;    /// requested alignment (at the moment only as check)
    PoolI!(ChunkGuard) poolChunks; /// pool for quick alloc of backing memory
    PoolI!(SegmentedArray!(T)) poolSegmented; /// pool for quick alloc of SegmentedArrays (without backing memory)
    KindRange kRange; /// range of the map (must be a subset of arrayStruct.kRange)
    size_t activeUsers=1;
    
    enum Flags{
        None,
        StandardBlock, // memory is a continuos block kindOffsets[0]..kindOffsets[$1], with each kind one after the
        // other without gaps
    }
    Flags flags=Flags.None;
    this(){}
    this(SegmentedArrayStruct arrayStruct,KindRange kRange=KindRange.all,size_t allocSize=0,size_t alignBytes=0,
        PoolI!(ChunkGuard) poolChunks=null, PoolI!(SegmentedArray!(T)) poolSegmented=null){
        this.arrayStruct=arrayStruct;
        arrayStruct.freeze();
        this.kRange=kRange.intersect(arrayStruct.kRange);
        this.allocSize=allocSize;
        this.alignBytes=alignBytes;
        this.poolChunks=poolChunks;
        this.poolSegmented=poolSegmented;
        auto nkinds=this.kRange.length;
        auto buf=new size_t[](2*nkinds);
        this.kindOffsets=buf[0..nkinds];
        auto kShift=this.kRange.kStart-arrayStruct.kRange.kStart;
        foreach(i,k;arrayStruct.kindStarts[kShift..kShift+nkinds]){
            this.kindOffsets[i]=(k-arrayStruct.kindStarts[kShift])*T.sizeof;
        }
        this.kindByteIncrements=buf[nkinds..2*nkinds];
        foreach(i,k;arrayStruct.kindIncrements[kShift..kShift+nkinds]){
            this.kindByteIncrements[i]=k*T.sizeof;
        }
        this.activeUsers=1;
        this.flags=Flags.StandardBlock;
        if (this.allocSize==0) {
            this.allocSize=(this.arrayStruct.kindStarts[kShift+nkinds]
                -this.arrayStruct.kindStarts[kShift])*T.sizeof;
        }
        if (this.poolChunks is null){
            if (this.allocSize>10*1024*1024){
                // use a global (non numa) pool
                this.poolChunks=new Pool!(ChunkGuard)(&allocBase);
            } else {
                this.poolChunks=cachedPool(&allocBase);
            }
        }
        if (this.poolSegmented is null){
            this.poolSegmented=cachedPool(&allocEmptySegArr);
        }
    }
    /// allocates a new block of storage (ChunkGuard)
    ChunkGuard allocBase(PoolI!(ChunkGuard)p){
        return new ChunkGuard(allocSize,alignBytes,typeHasPointers!(T)(),p);
    }
    /// allocates a new segmented array (without backing memory)
    SegmentedArray!(T) allocEmptySegArr(PoolI!(SegmentedArray!(T))p){
        return new SegmentedArray!(T)(p);
    }
    /// new array (with data)
    SegmentedArray!(T) newArray(){
        auto guard=poolChunks.getObj();
        auto sArr=poolSegmented.getObj();
        sArr.reset(this,guard);
        guard.release();
        return sArr;
    }
    /// new array (with a subset of the data)
    SegmentedArray!(T) newArray(KindRange kr){
        // realloc another map is it is a much smaller subset???
        auto guard=poolChunks.getObj();
        auto sArr=poolSegmented.getObj();
        sArr.reset(this,guard,kr);
        guard.release();
        return sArr;
    }
    /// new array (with a BulkArray)
    SegmentedArray!(T) newArray(BulkArray!(T) bArr,KindRange kr=KindRange.all){
        auto sArr=poolSegmented.getObj();
        sArr.reset(this,bArr,kr);
        return sArr;
    }
    /// if the array is contiguous
    bool contiguous(){
        return (flags&Flags.StandardBlock)!=0;
    }
    // fully disallow serialization, as it is useful mostly for debugging purposes...
    mixin(serializeSome("dchem.SegArrMemMap!("~T.stringof~")",`Describes the memory mapping of a SegmentedArray`,
        `arrayStruct
        kindOffsets
        kindByteIncrements
        allocSize
        alignBytes`));
    mixin printOut!();
    
    // pool methods
    
    /// returns a new segmented array
    SegmentedArray!(T) getObj(){
        return newArray();
    }
    /// returns an instance, so that it can be reused
    void giveBack(SegmentedArray!(T) obj){
        if (obj is null) return;
        if (obj.guard !is null){
            obj.guard.release;
            obj.guard=null;
        }
        if (obj.pool!is null){
            obj.pool.giveBack(obj);
        } else {
            // actively delete??
        }
    }
    /// should discard all the cached objects (no guarantee)
    void flush(){
        poolChunks.flush();
        poolSegmented.flush();
    }
    /// should not cache any objects from now on (no guarantee)
    /// this is suboptimal if the pools are shared, use some reference counting or avoid propagating???
    void stopCaching(){
        poolChunks.stopCaching();
        poolSegmented.stopCaching();
    }
    /// add an active user (when created one active user is automatically added)
    void addUser(){
        if (atomicAdd(activeUsers,cast(size_t)1)==0){
            throw new Exception("addUser called on non used pool",__FILE__,__LINE__);
        }
    }
    /// removes an active user (if there are 0 active users stopCaching is called)
    void rmUser(){
        auto oldUsers=atomicAdd(activeUsers,-cast(size_t)1);
        if (oldUsers==0){
            throw new Exception("rmUser called on non used pool",__FILE__,__LINE__);
        }
        if (oldUsers==1){
            stopCaching();
        }
    }
    /// tries to reuse the pools from the given SegArrMemMap
    void consolidate(V)(SegArrMemMap!(V) baseVal){
        if (baseVal is null) return;
        static if (is(V==T)){
            baseVal.poolSegmented.addUser();
            poolSegmented.rmUser();
            poolSegmented=baseVal.poolSegmented;
        }
        if (baseVal.allocSize==allocSize && baseVal.alignBytes>=alignBytes){
            baseVal.poolChunks.addUser();
            poolChunks.rmUser();
            poolChunks=baseVal.poolChunks;
        }
    }
}

/// segmented (level/kinds) array (with kind uniform dimension)
final class SegmentedArray(T){
    void *basePtr; /// this is guard.dataPtr
    size_t[] kindOffsets; /// offsets wrt. to basePtr (in bytes)
    size_t[] kindByteIncrements; /// increments for each particle (in bytes)
    index_type[] mKindDims; /// number of basic elements in each particle (at least 1 if Min1)
    SegArrMemMap!(T) arrayMap;
    ChunkGuard guard;
    KindRange  kRange;
    alias T dtype;
    alias BulkArray!(T).basicDtype basicDtype;
    static size_t defaultOptimalBlockSize=32*1024/T.sizeof;
    PoolI!(SegmentedArray) pool;
    
    mixin(serializeSome("dchem.sys.SegmentedArray!("~T.stringof~")","A segmented array i.e some particle associated property.",
        "kRange|kindOffsets|kindByteIncrements|mKindDims|arrayMap|guard"));
    mixin printOut!();
    
    /+BulkArray!(T) data(){
        return _data[kindStarts[0],kindStarts[$-1]];
    }+/

    SegmentedArrayStruct arrayStruct(){
        return arrayMap.arrayStruct;
    }
    /// initializes with map, guard & range
    void reset(SegArrMemMap!(T) arrayMap,ChunkGuard guard,KindRange kr){
        assert(this.arrayMap is null && this.guard is null && this.basePtr is null,"reset called on non cleared SegmentedArray");
        this.arrayMap=arrayMap;
        this.guard=guard.retain(); // allow null guard???
        this.basePtr=guard.dataPtr;
        assert(kr in arrayMap.kRange,"kr out of range of map");
        this.kRange=kr;
        auto kStart=arrayMap.kRange.kStart;
        this.kindOffsets=arrayMap.kindOffsets[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
        this.kindByteIncrements=arrayMap.kindByteIncrements[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
        this.mKindDims=arrayMap.arrayStruct.mKindDims[(this.kRange.kStart-arrayMap.arrayStruct.kRange.kStart)..
            (this.kRange.kEnd-arrayMap.arrayStruct.kRange.kStart)]; // duplicate this in the map?
    }
    /// initializes with map & guard
    void reset(SegArrMemMap!(T) arrayMap,ChunkGuard guard){
        reset(arrayMap,guard,arrayMap.kRange);
    }
    /// initializes with a BulkArray, that must be compatible with the memory layout of arrayMap
    void reset(SegArrMemMap!(T) arrayMap,BulkArray!(T)bArr,KindRange kr=KindRange.all){
        assert(this.arrayMap is null && this.guard is null && this.basePtr is null,"reset called on non cleared SegmentedArray");
        this.arrayMap=arrayMap;
        guard=bArr.guard;
        if (guard is null){
            // copy...
            guard=arrayMap.poolChunks.getObj();
            this.basePtr=guard.dataPtr;
            this.kRange=kr.intersect(arrayMap.kRange);
            auto kStart=arrayMap.kRange.kStart;
            auto kStart2=arrayMap.arrayStruct.kRange.kStart;
            this.kindOffsets=arrayMap.kindOffsets[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
            this.kindByteIncrements=arrayMap.kindByteIncrements[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
            this.mKindDims=arrayMap.arrayStruct.mKindDims[(this.kRange.kStart-kStart2)..(this.kRange.kEnd-kStart2)];
            assert(bArr.length==dataLength);
            foreach(i,ref v;sDataLoop){ // could be parallel...
                v=bArr[i];
            }
        } else {
            // check compatibility
            guard.retain();
            this.basePtr=guard.dataPtr;
            this.kRange=kr.intersect(arrayMap.kRange);
            auto kStart=arrayMap.kRange.kStart;
            auto kStart2=arrayMap.arrayStruct.kRange.kStart;
            this.kindOffsets=arrayMap.kindOffsets[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
            this.kindByteIncrements=arrayMap.kindByteIncrements[(this.kRange.kStart-kStart)..(this.kRange.kEnd-kStart)];
            this.mKindDims=arrayMap.arrayStruct.mKindDims[(this.kRange.kStart-kStart2)..(this.kRange.kEnd-kStart2)];
            assert(bArr.length==dataLength);
            auto bArr2=support();
            if (bArr2.ptr!is bArr.ptr || bArr2.length!=bArr.length){
                throw new Exception("bArr incompatible with arrayMap",__FILE__,__LINE__);
            }
        }
    }
    

    // internal for serialization
    this(){ }
    this(PoolI!(SegmentedArray) pool){
        pool=pool;
    }
    this(SegArrMemMap!(T) arrayMap,ChunkGuard guard,KindRange kr,PoolI!(SegmentedArray) pool=null){
        this(pool);
        reset(arrayMap,guard,kr);
    }
    this(SegArrMemMap!(T) arrayMap,ChunkGuard guard,PoolI!(SegmentedArray) pool=null){
        this(pool);
        reset(arrayMap,guard);
    }
    this(SegArrMemMap!(T) arrayMap,BulkArray!(T)bArr,PoolI!(SegmentedArray) pool=null){
        this(pool);
        reset(arrayMap,bArr);
    }
    /// if the array is stored contiguosly
    bool contiguous(){
        return arrayMap.contiguous;
    }
    void clear(){
        if (guard!is null){
            guard.release;
            guard=null;
            basePtr=null;
        }
        arrayMap=null;
        kindOffsets=[];
        kindByteIncrements=[];
        mKindDims=[];
    }
    /// returns a bulk array that can support the current range
    BulkArray!(T)support(){
        if (kindOffsets.length==0) return BulkArray!(T)();
        auto kEnd=kRange.kEnd-kRange.kStart;
        auto aStruct=arrayMap.arrayStruct;
        auto kShift=kRange.kStart-aStruct.kRange.kStart;
        size_t minS=kindOffsets[0];
        size_t maxS=kindOffsets[0]+(aStruct.kindStarts[kShift+1]-aStruct.kindStarts[kShift])*T.sizeof;
        for (auto i=1;i<kEnd;++i){
            if (maxS<=kindOffsets[i]){ // assumes non overlapping segments
                maxS=kindOffsets[i]+(aStruct.kindStarts[i+kShift+1]-aStruct.kindStarts[i+kShift])*T.sizeof;
            }
            if (minS>kindOffsets[i]) minS=kindOffsets[i];
        }
        assert(((maxS-minS)%T.sizeof)==0,"non aligned offsets");
        return BulkArray!(T)((cast(T*)(basePtr+minS))[0..((maxS-minS)/T.sizeof)],guard);
    }
/+    /// allocates a new SegmentedArray with the given kind dimensions
    /// min1 
    this(SegmentedArrayStruct arrayStruct, BulkArray!(T)data=BulkArray!(T).dummy,
        KindRange kRange=KindRange.all,index_type[] kindStarts=null,PoolI!(SegmentedArray) pool=null)
    in{
        assert(arrayStruct!is null,"arrayStruct must be valid");
        auto myKRange=kRange;
        if (kRange.kEnd==KindIdx.init) myKRange=arrayStruct.kRange;
        if (kindStarts!is null){
            arrayStruct.freeze;
            auto nkinds=cast(size_t)(myKRange.kEnd-myKRange.kStart);
            assert(myKRange in arrayStruct.kRange,"kRange out of bounds");
            assert(kindStarts.length==nkinds+1,"kindStarts has wrong size");
            auto kindShift=cast(size_t)(myKRange.kStart-arrayStruct.kRange.kStart);
            for (size_t i=0;i<=nkinds;++i){
                assert(kindStarts[i]==kindStarts[0]+arrayStruct.kindStarts[kindShift+i]-arrayStruct.kindStarts[kindShift],
                    "invalid kindStarts");
            }
        }
    } body {
        assert(arrayStruct !is null,"arrayStruct must be valid");
        this.arrayStruct=arrayStruct.freeze;
        this.kRange=kRange;
        if (kRange.kEnd==KindIdx.init){
            assert(kRange.kStart==0);
            this.kRange=arrayStruct.kRange;
        }
        assert(this.kRange in arrayStruct.kRange);
        this.kindStarts=kindStarts;
        if (this.kindStarts.length==0){
            auto nkinds=cast(size_t)(this.kRange.kEnd-this.kRange.kStart);
            this.kindStarts=new index_type[nkinds+1];
            auto kindShift=cast(size_t)(this.kRange.kStart-arrayStruct.kRange.kStart);
            for (size_t i=0;i<=nkinds;++i){
                this.kindStarts[i]=arrayStruct.kindStarts[kindShift+i]-arrayStruct.kindStarts[kindShift];
            }
        }
        _data=data;
        if (BulkArrayIsDummy(data)){
            _data=BulkArray!(T)(this.kindStarts[this.kindStarts.length-1]-this.kindStarts[0]);
        }
        this.pool=pool;
    }+/
    /// array of elements in the given range
    SegmentedArray!(T) opIndex(KindRange kr){
        auto res=arrayMap.poolSegmented.getObj();
        KindRange krCommon=kRange.intersect(kr);
        res.reset(arrayMap,guard,krCommon);
        return res;
    }
    /// array of elements for local kind k
    BulkArray!(T) opIndex(KindIdx k){
        assert(k in kRange,"kind out of range");
        auto aStruct=arrayMap.arrayStruct;
        auto res=BulkArray!(T)((cast(T*)(basePtr+kindOffsets[k-kRange.kStart]))[0..(aStruct.kindStarts[k-aStruct.kRange.kStart+1]-
            aStruct.kindStarts[k-aStruct.kRange.kStart])],guard);
        assert(res.ptr>=basePtr && res.ptr+res.length<=(basePtr+guard.dataLen),"data out of range");
        return res;
    }
    /// gets the particle using the local particle numbering, but that might be outside the kind range
    /// could be optimized more, but will never be very efficient...
    T[] getMaybeInRange(LocalPIndex p){
        auto kindIdx=p.kind();
        if (kindIdx in kRange){
            auto kStart=arrayMap.arrayStruct.submapping.lKRange.kStart;
            if (p.particle() < arrayMap.arrayStruct.nParticles(kindIdx)){
                auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+p.particle()*kindByteIncrements[kindIdx-kRange.kStart]);
                auto res=startPtr[0..mKindDims[kindIdx-kRange.kStart]];
                assert(res.ptr>=basePtr && res.ptr+res.length<(basePtr+guard.dataLen),"data out of range");
            }
        }
        return null;
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    T[] opIndex(LocalPIndex p){
        auto kindIdx=p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
        assert(p.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range"); // make this test quicker?
        auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+p.particle()*kindByteIncrements[kindIdx-kRange.kStart]);
        auto res=startPtr[0..mKindDims[kindIdx-kRange.kStart]];
        assert(res.ptr>=basePtr && res.ptr+res.length<=(basePtr+guard.dataLen),"data out of range");
        return res;
    }
    /// array of elements for a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T[] opIndex(PIndex p){
        if ((arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0){
            auto kindIdx=p.kind();
            auto pos=p.particle();
            assert(kindIdx in kRange,"kind out of range");
            assert(pos < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range");
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+pos*kindByteIncrements[kindIdx-kRange.kStart]);
            auto res=startPtr[0..mKindDims[kindIdx-kRange.kStart]];
            assert(res.ptr>=basePtr && res.ptr+res.length<=(basePtr+guard.dataLen),"data out of range");
            return res;
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            assert(l.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range, internal error");
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+l.particle()*kindByteIncrements[kindIdx-kRange.kStart]);
            auto res=startPtr[0..mKindDims[kindIdx-kRange.kStart]];
            assert(res.ptr>=basePtr && res.ptr+res.length<=(basePtr+guard.dataLen),"data out of range");
            return res;
        }
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T *ptrI(LocalPIndex p,index_type i){
        auto kindIdx=p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
        assert(i>=0 && (i<arrayStruct.mKindDims[kindIdx-kRange.kStart]||
                (i==0&&(arrayStruct.flags & SegmentedArrayStruct.Flags.Min1)!=0)),"index i out of bounds");
        assert(p.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range, internal error");
        assert(i<mKindDims[kindIdx-kRange.kStart],"i out of range"); // always allow i=0??
        auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+
            p.particle()*kindByteIncrements[kindIdx-kRange.kStart]+i*T.sizeof);
        return startPtr;
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T opIndex(LocalPIndex p,index_type i){
        return *ptrI(p,i);
    }
    /// sets the value for element i of the particle using the local particle numbering
    /// (has to be in the kind range)
    void opIndexAssign(T val,LocalPIndex p,index_type i){
        static if (is(typeof(delegate void(){ val[]=val; }))){
            (*ptrI(p,i))[]=val;
        } else {
            *ptrI(p,i)=val;
        }
    }
    /// array of elements for a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    void opIndexAssign(T val,PIndex p,index_type i){
        if ((arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0){
            auto kindIdx=p.kind();
            auto pos=p.particle();
            assert(kindIdx in kRange,"kind out of range");
            assert(pos < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range");
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+pos*kindByteIncrements[kindIdx-kRange.kStart]);
            assert(i>=0&&i<mKindDims[kindIdx-kRange.kStart],"in particle index i out of bounds");
            startPtr[i]=val;
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            assert(l.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range, internal error");
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+l.particle()*kindByteIncrements[kindIdx-kRange.kStart]);
            assert(i>=0&&i<mKindDims[kindIdx-kRange.kStart],"in particle index i out of bounds");
            startPtr[i]=val;
        }
    }
    /// address of element i of a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T *ptrI(PIndex p,index_type i){
        if ((arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0){
            auto kindIdx=p.kind();
            assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
            assert(p.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range, internal error");
            assert(i>=0 && i<mKindDims[kindIdx-kRange.kStart],"i out of range"); // always allow i=0??
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+
                p.particle()*kindByteIncrements[kindIdx-kRange.kStart]+i*T.sizeof);
            return startPtr;
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
            assert(l.particle() < arrayMap.arrayStruct.nParticles(kindIdx),"particle index out of range, internal error");
            assert(i>=0 && i<mKindDims[kindIdx-kRange.kStart],"i out of range"); // always allow i=0??
            auto startPtr=cast(T*)(basePtr+kindOffsets[kindIdx-kRange.kStart]+
                l.particle()*kindByteIncrements[kindIdx-kRange.kStart]+i*T.sizeof);
            return startPtr;
        }
    }
    /// element i of a (global) particle p array
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T opIndex(PIndex p,index_type i){
        return *ptrI(p,i);
    }
    void opSliceAssign(V)(V v){
        static if (is(typeof(opSliceAssignEl(v)))){
            opSliceAssignEl(v);
        } else static if (is(typeof(this.opSliceAssignT(v)))){
            this.opSliceAssignT(v);
        } else static if (is(typeof(v.copyTo(this)))){
            v.copyTo(this);
        } else {
            static assert(0,"cannot assign from "~V.stringof~" to SegmentedArray!("~T.stringof~")");
        }
    }
    /// copies from an array to this
    void opSliceAssignT(V)(SegmentedArray!(V) val){
        static if (is(T==V)) {
            if (val is this) return;
	}
        assert(arrayStruct==val.arrayStruct,"different structs");
        assert(kRange==val.kRange,"different kRanges");
        int kEnd=kRange.kEnd-kRange.kStart;
        auto kStarts=arrayMap.arrayStruct.kindStarts[kRange.kStart-arrayMap.arrayStruct.kRange.kStart..$];
        for (int ik=0;ik<kEnd;++ik){
            auto startPtr1=cast(T*)(basePtr+kindOffsets[ik]);
            auto startPtr2=cast(V*)(val.basePtr+val.kindOffsets[ik]);
            static if (is(T==V)){
		assert(kindByteIncrements[ik]==val.kindByteIncrements[ik],"different particle sizes");
                memcpy(startPtr1,startPtr2,T.sizeof*(kStarts[ik+1]-kStarts[ik]));
            } else {
                for(size_t p=kStarts[ik+1]-kStarts[ik];p!=0;--p){
                    static if(is(typeof((*startPtr1)[]=(*startPtr2)))){
                        (*startPtr1)[]=(*startPtr2);
                    } else static if(is(V:T)){
                        *startPtr1=cast(T)(*startPtr2);
                    } else {
                        *startPtr1=convertTo!(T)(*startPtr2);
                    }
                    ++startPtr1;
                    ++startPtr2;
                }
            }
        }
    }
    /// copies from an array to this
    void opSliceAssignEl(T val){
        int kEnd=kRange.kEnd-kRange.kStart;
        auto kStarts=arrayMap.arrayStruct.kindStarts[kRange.kStart-arrayMap.arrayStruct.kRange.kStart..$];
        for (int ik=0;ik<kEnd;++ik){
            auto startPtr1=cast(T*)(basePtr+kindOffsets[ik]);
            for(size_t p=kStarts[ik+1]-kStarts[ik];p!=0;--p){
                static if(is(typeof((*startPtr1)[]=val))){
                    (*startPtr1)[]=val;
                } else {
                    *startPtr1=val;
                }
                ++startPtr1;
            }
        }
    }
    /// copies this array to the given SegmentedArray, tryig to reuse its memory allocations
    void dupTo(V)(SegmentedArray!(V) val){
        static if(is(T==V)){
            if (val.kRange!=kRange ||
		val.arrayMap.arrayStruct != arrayMap.arrayStruct){
		val.clear();
                // avoid using this mapping if it is a much smaller subset?
                val.reset(arrayMap,arrayMap.poolChunks.getObj());
                if (val.guard!is null) val.guard.release();
            }
            val.opSliceAssignT!(T)(this);
        } else {
            if (val.kRange!=kRange||
		val.arrayMap.arrayStruct != arrayMap.arrayStruct){
                val.clear();
                auto newMap=new SegArrMemMap!(V)(arrayStruct);
                newMap.consolidate(arrayMap);
                val.reset(newMap,newMap.poolChunks.getObj());
                if (val.guard!is null) val.guard.release();
            }
            val.opSliceAssignT!(T)(this);
        }
    }
    /// returns a copy of the segmented array
    SegmentedArray!(V) dupT(V=T)(){
        static if (is(V==T)){
            auto res=arrayMap.newArray();
            dupTo(res);
            return res;
        } else {
            auto newMap=new SegArrMemMap!(V)(arrayStruct);
            newMap.consolidate(arrayMap);
            auto res=newMap.newArray();
            this.dupTo(res);
            return res;
        }
    }
    SegmentedArray dup(){
        return dupT!()();
    }
    /// returns a copy of the segmented array
    SegmentedArray deepdup(){
        assert(0,"unimplemented");
        return null;
    }
    
    /// length of a particle based loop on this (i.e skipping kinds without particles, and counting
    /// the number of particles for dimension of size 0 if Min1, i.e. a single value per kind)
    /// this is the size of sLoop and pLoop. The real storage size is .dataLength
    size_t length(){
        size_t iterRef=this.dataLength;
        auto submap=this.arrayStruct.submapping;
        if((arrayStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0){
            for (auto k=this.kRange.kStart;k<this.kRange.kEnd;++k){
                if (this.arrayStruct.kindDim(k)==0){
                    iterRef-=1;
                    iterRef+=submap.kindStarts[k-submap.lKRange.kStart+1]-submap.kindStarts[k-submap.lKRange.kStart];
                }
            }
        }
        return iterRef;
    }
    /// returns the length of the data actually stored (this is the "dimension" of this vector in T units)
    size_t dataLength(){
        auto aStruct=arrayMap.arrayStruct;
        auto kStart=aStruct.kRange.kStart;
        auto res=aStruct.kindStarts[kRange.kEnd-kStart]-aStruct.kindStarts[kRange.kStart-kStart];
        return res;
    }
    
    /// loops on the particles and each element of a multivaluated segmented array.
    /// for consistency always skips kinds without particles (even if Min1)
    struct PLoop(int pFlags){
        size_t optimalBlockSize;
        SegmentedArray array;
        int opApply(int delegate(ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","if (auto res=dlg(*arrayPtr)){ *finalRes=res; return; }","",
            
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayMKindDims;++ii){
                if (auto res=dlg(arrayPtr[ii])){ *finalRes=res; return; }
            }`,""]));
            return result;
        }
        static if ((pFlags&ParaFlags.DataLoop)!=0){
            int opApply(int delegate(ref size_t gIdx,ref T)loopBody){
                int result=0;
                mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
                "int delegate(ref size_t,ref T) dlg; int* finalRes; SegmentedArrayStruct aStruct;","",
                "mainContext.dlg=loopBody;mainContext.finalRes=&result;mainContext.aStruct=array.arrayStruct;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
                ["if ((*finalRes)!=0) return; size_t gIdx=aStruct.kindStarts[kind-aStruct.kRange.kStart]+start;","if (auto res=dlg(gIdx,*arrayPtr)){ *finalRes=res; return; } ++gIdx;","",

                "if ((*finalRes)!=0) return; size_t gIdx=aStruct.kindStarts[kind-aStruct.kRange.kStart]+start;",`
                for (size_t ii=0;ii<arrayMKindDims;++ii){
                    if (auto res=dlg(gIdx,arrayPtr[ii])){ *finalRes=res; return; }
                    ++gIdx;
                }`,""]));
                return result;
            }
        } else {
            int opApply(int delegate(ref LocalPIndex lIdx,ref T[])loopBody){
                int result=0;
                mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
                "int delegate(ref LocalPIndex,ref T[]) dlg; int* finalRes;","",
                "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
                ["if ((*finalRes)!=0) return;","auto myArr=arrayPtr[0..1];if (auto res=dlg(localPIndex,myArr)){ *finalRes=res; return; }","",
            
                "if ((*finalRes)!=0) return;",`
                auto myArr=arrayPtr[0..arrayMKindDims];
                if (auto res=dlg(localPIndex,myArr)){ *finalRes=res; return; }
                /+for (size_t ii=0;ii<arrayMKindDims;++ii){
                    if (auto res=dlg(localPIndex,arrayPtr[ii])){ *finalRes=res; return; }
                }+/`,""]));
                return result;
            }
        }
        int opApply(int delegate(ref PIndex pIdx,ref LocalPIndex lIdx,ref T[])loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref PIndex, ref LocalPIndex, ref T[]) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","auto myArr=arrayPtr[0..1]; if (auto returnV=dlg(*pIndexPtr,localPIndex,myArr)){ *finalRes=returnV; return; }","",
            "if ((*finalRes)!=0) return;",`
            auto myArr=arrayPtr[0..arrayMKindDims];
            if (auto returnV=dlg(*pIndexPtr,localPIndex,myArr)){ *finalRes=returnV; return; }
            /+for (size_t ii=0;ii<arrayMKindDims;++ii){
                if (auto returnV=dlg(*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
            }+/`,""]));
           return result;
        }
        int opApply(int delegate(ref size_t i,ref PIndex pIdx,ref LocalPIndex lIdx,ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref size_t i,ref PIndex, ref LocalPIndex, ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","size_t ii=0; if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,*arrayPtr)){ *finalRes=returnV; return; }","",
            
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayMKindDims;++ii){
                if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
            }`,""]));
           return result;
        }
    }
    /// full parallel loop
    PLoop!(ParaFlags.FullPara) pLoop(size_t optSize=defaultOptimalBlockSize){
        PLoop!(ParaFlags.FullPara) res;
        res.optimalBlockSize=optSize;
        res.array=this;
        return res;
    }
    /// loop parallel only between kinds
    PLoop!(ParaFlags.KindPara) pLoopKinds(size_t optSize=defaultOptimalBlockSize){
        PLoop!(ParaFlags.KindPara) res;
        res.optimalBlockSize=optSize;
        res.array=this;
        return res;
    }
    /// sequential loop
    PLoop!(ParaFlags.Sequential) sLoop(){
        PLoop!(ParaFlags.Sequential) res;
        res.optimalBlockSize=defaultOptimalBlockSize;
        res.array=this;
        return res;
    }
    /// parallel data loop
    PLoop!(ParaFlags.FullPara|ParaFlags.DataLoop) pDataLoop(size_t optSize=defaultOptimalBlockSize){
        PLoop!(ParaFlags.FullPara|ParaFlags.DataLoop) res;
        res.optimalBlockSize=optSize;
        res.array=this;
        return res;
    }
    /// sequential data loop
    PLoop!(ParaFlags.Sequential|ParaFlags.DataLoop) sDataLoop(){
        PLoop!(ParaFlags.Sequential|ParaFlags.DataLoop) res;
        res.optimalBlockSize=defaultOptimalBlockSize;
        res.array=this;
        return res;
    }

    /// opBypax
    void opBypax(V)(SegmentedArray!(V) x,basicDtype a=1,basicDtype b=1){
        if (x.arrayStruct is arrayStruct && x.arrayMap.contiguous && arrayMap.contiguous){
            scope a1=a2NA(support.basicData);
            scope a2=a2NA(x.support.basicData);
            a1.opBypax(a2,a,b);
            return;
        }
        auto optimalBlockSize=defaultOptimalBlockSize;
        auto y=this;
        if(b==1){
            if(a==1){
                mixin(segArrayMonoLoop(ParaFlags.FullPara/+|ParaFlags.DataLoop+/,"iterContext",["y","x"],"","",
                "",`
                visitKind=visitKind&&(!outOfRange);
                assert((newK.xMKindDims<=1&&newK.yMKindDims<=1)||newK.xMKindDims==newK.yMKindDims,"variable combination of MKindDims not implemented in opBypax");`,"",
                ["","*yPtr += convertTo!(T)(*xPtr);","",
                
                "",`
                for (size_t ii=0;ii<xMKindDims;++ii){
                    yPtr[ii] += convertTo!(T)(xPtr[ii]);
                }`,""]));
            } else {
                alias optimalBlockSize optimalBlockSize_1;
                mixin(segArrayMonoLoop(ParaFlags.FullPara/+|ParaFlags.DataLoop+/,"iterContext2",["y","x"],
                "basicDtype a;","_1",
                "mainContext_1.a=a;",`
                visitKind_1=visitKind_1&&(!outOfRange_1);
                assert((newK_1.xMKindDims<=1&&newK_1.yMKindDims<=1)||newK_1.xMKindDims==newK_1.yMKindDims,"variable combination of MKindDims not implemented in opBypax");`,"",
                ["","*yPtr += convertTo!(T)((*xPtr)*this.a);","",
                
                "",`
                for (size_t ii=0;ii<xMKindDims;++ii){
                    yPtr[ii] += convertTo!(T)(xPtr[ii]*this.a);
                }`,""]));
            }
        } else if(b==0){
            alias optimalBlockSize optimalBlockSize_2;
            mixin(segArrayMonoLoop(ParaFlags.FullPara/+|ParaFlags.DataLoop+/,"iterContext_2",["y","x"],
            "basicDtype a;","_2",
            "mainContext_2.a=a;",`
            visitKind_2=visitKind_2&&(!outOfRange_2);
            assert((newK_2.xMKindDims<=1&&newK_2.yMKindDims<=1)||newK_2.xMKindDims==newK_2.yMKindDims,"variable combination of MKindDims not implemented in opBypax");`,"",
            ["","*yPtr = convertTo!(T)((*xPtr)*a);","",
            
            "",`
            for (size_t ii=0;ii<xMKindDims;++ii){
                yPtr[ii] = convertTo!(T)(xPtr[ii]*a);
            }`,""]));
        } else {
            alias optimalBlockSize optimalBlockSize_3;
            mixin(segArrayMonoLoop(ParaFlags.FullPara/+|ParaFlags.DataLoop+/,"iterContext_3",["y","x"],
            "basicDtype a;basicDtype b;","_3",
            "mainContext_3.a=a;mainContext_3.b=b;",`
            visitKind_3=visitKind_3&&(!outOfRange_3);
            assert((newK_3.xMKindDims<=1&&newK_3.yMKindDims<=1)||newK_3.xMKindDims==newK_3.yMKindDims,"variable combination of MKindDims not implemented in opBypax");`,"",
            ["","*yPtr = convertTo!(T)((*yPtr)*b+(*xPtr)*a);","",
            
            "",`
            for (size_t ii=0;ii<xMKindDims;++ii){
                yPtr[ii] = convertTo!(T)(yPtr[ii]*b+xPtr[ii]*a);
            }`,""]));
        }
    }
    
    static if (is(typeof(basicDtype.init*basicDtype.init))){
        void opMulAssign()(basicDtype scale){
            if (contiguous){
                scope a=a2NA(this.support.basicData);
                a*=scale;
            } else {
                auto optimalBlockSize=defaultOptimalBlockSize;
                auto x=this;
                mixin(segArrayMonoLoop(ParaFlags.FullPara|ParaFlags.DataLoop,"iterContext",["x"],
                "basicDtype a;","",
                "mainContext.a=scale;",`
                visitKind=visitKind&&(!outOfRange);`,"",
                ["","for(size_t mIdx=0;mIdx<T.sizeof/basicDtype.sizeof;++mIdx){ (cast(basicDtype*)xPtr)[mIdx] *= a; }","",

                "",`
                for (size_t ii=0;ii<xMKindDims*(T.sizeof/basicDtype.sizeof);++ii){
                    (cast(basicDtype*)xPtr)[ii] *= a;
                }`,""]));
            }
        }
    }
    void opMulAssign(V)(SegmentedArray!(V) y){
        if (x.arrayStruct is arrayStruct){
            scope a1=a2NA(data.basicData);
            scope a2=a2NA(y.data.basicData);
            a1*=a2;
            return;
        }
        auto optimalBlockSize=defaultOptimalBlockSize;
        auto x=this;
        static if (is(typeof(function(T t){ t *= t; }))){
            mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext_2",["x","y"],
            "T scale;","",
            "mainContext.scale=scale;",`
            visitKind=visitKind&&(!outOfRange);
            assert((newK.xMKindDims<=1&&newK.yMKindDims<=1)||newK.xMKindDims==newK.yMKindDims,"variable combination of MKindDims not implemented in opMulAssign");`,"",
            ["","*xPtr *= *yPtr;","",
            
            "",`
            for (size_t ii=0;ii<xMKindDims;++ii){
                xPtr[ii] *= convertTo!(T)(yPtr[ii]);
            }`,""]));
        } else {
            alias optimalBlockSize optimalBlockSize_2;
            mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext_2",["x","y"],
            "T scale;","_2",
            "mainContext_2.scale=scale;",`
            visitKind_2=visitKind_2&&(!outOfRange_2);
            assert((newK_2.xMKindDims<=1&&newK_2.yMKindDims<=1)||newK_2.xMKindDims==newK_2.yMKindDims,"variable combination of MKindDims not implemented in opMulAssign");`,"",
            ["","*xPtr = (*xPtr)*(*yPtr);","",
            
            "",`
            for (size_t ii=0;ii<xMKindDims;++ii){
                xPtr[ii] = convertTo!(T)(xPtr[ii]*yPtr[ii]);
            }`,""]));
        }
    }
    /// an array just like this, but with unitialized values
    SegmentedArray emptyCopy(){
        return arrayMap.newArray(kRange);
    }
    /// gives back the current array, its use after this call is an error
    void giveBack(){
        if (guard!is null){
            guard.release();
            guard=null;
        }
        if (pool!is null){
            pool.giveBack(this);
        }
    }
}

/// inner loop variant (exec* method)
char[] segArrayContextExecStr(int pFlags,char[] iterName, char[][]namesLocal,char[] contextExtra,char[] uniq,char[] execN,
    char[]pVisitLocalStart,char[]pVisit,char[]pVisitLocalEnd){
    char[]visitorStr=`
        void exec`~execN~`(){
            if (this.context.exception !is null) return;`;
    visitorStr~="\n";
    visitorStr~=pVisitLocalStart;
    visitorStr~="\n";
    if (pFlags==ParaFlags.FullPara){
        visitorStr~=`
            if (this.end-this.start>this.optimalBlockSize*3/2){
                auto mid=cast(ParticleIdx)((this.end-this.start)/2);
                if (mid>this.optimalBlockSize){ // tries to make optimalBlockSize a possible fast path
                    mid=cast(ParticleIdx)(((mid+this.optimalBlockSize-1)/this.optimalBlockSize)*this.optimalBlockSize);
                }
                auto firstHalf=this.alloc();
                firstHalf.end=this.start+mid;
                Task("SegArrLoop1`~iterName~`",&firstHalf.exec`~execN~`).appendOnFinish(&firstHalf.giveBack).autorelease.submit; // would submitYield be better (less suspended tasks, but more suspensions)??? 
                auto secondHalf=this.alloc();
                secondHalf.start=this.start+mid;`;
        foreach (name;namesLocal){
            visitorStr~=`
                secondHalf.`~name~`PtrStart=cast(typeof(secondHalf.`~name~`PtrStart))(
                    cast(size_t)secondHalf.`~name~`PtrStart+this.`~name~`ByteIncrements*mid);`;
        }
        visitorStr~=`
                Task("SegArrLoop2`~iterName~`",&secondHalf.exec`~execN~`).appendOnFinish(&secondHalf.giveBack).autorelease.submit;
            } else {`;
    } else {
        visitorStr~=`
            {`;
    }
    visitorStr~=`
                try{
                    LocalPIndex localPIndex=LocalPIndex(this.kind,this.start);
                    PIndex * pIndexPtr=this.pIndexPtrStart;`;
    foreach (name;namesLocal){
        visitorStr~=`
                    `~name~`.dtype* `~name~`Ptr=this.`~name~`PtrStart;`;
    }
    visitorStr~=`
                    for (ParticleIdx index=this.start;index<this.end;++index){
                        this.lIndex=index;
                        `~pVisit~`
                        ++localPIndex;
                        ++pIndexPtr;`;
    foreach (name;namesLocal){
        visitorStr~=`
                        `~name~`Ptr=cast(typeof(`~name~`Ptr))(
                            cast(size_t)`~name~`Ptr+`~name~`ByteIncrements);`;
    }
    visitorStr~=`
                    }
                } catch (Exception e){
                    this.context.exception=e;
                }
            }
            `~pVisitLocalEnd~`
        }`;
    return visitorStr;
}

/// the whole context struct, can have several local visitors (pVisitors), those are
/// always sets of 3 strings: pVisitLocalStart,pVisit,pVisitLocalEnd, and they create
/// each one exec context: exec0,exec1,...
char[] segArrayContextStr(int pFlags,char[] iterName, char[][]namesLocal,char[] contextExtra,
    char[] uniq,char[][] pVisitors=[])
{
    assert(pVisitors.length%3==0,"otherVisitors is supposed to be a multiple of 3");
    char[] visitorStr=`
    struct `~iterName~`{`;
    foreach (name;namesLocal){
        visitorStr~=`
        `~name~`.dtype* `~name~`PtrStart;
        size_t `~name~`MKindDims;
        size_t `~name~`ByteIncrements;`;
    }
    visitorStr~=`
        size_t optimalBlockSize;
        KindIdx kind;
        ParticleIdx start;
        ParticleIdx lIndex;
        ParticleIdx end;
        PIndex *pIndexPtrStart;
        size_t maxInnerDim;
        Exception exception;
        `~iterName~`* context;
        `~iterName~`* next;`;
    visitorStr~="\n";
    visitorStr~=contextExtra;
    visitorStr~=`
        typeof(this) alloc(){
            auto res=popFrom(this.context.next);
            if (res is null){
                res=new `~iterName~`;
            }
            *res=*this;
            res.next=this;
            return res;
        }`;
    for(int i=0;i<pVisitors.length;i+=3){
        visitorStr~=segArrayContextExecStr(pFlags,iterName,namesLocal,contextExtra,uniq,ctfe_i2a(i/3),
            pVisitors[i],pVisitors[i+1],pVisitors[i+2]);
    }
    visitorStr~=`
        void giveBack(){
            insertAt(this.context.next,this);
        }
    }
    `;
    return visitorStr;
}

/// loop on the kinds (main external loop)
char[] segArrayKLoopStr(int pFlags,char[] iterName, char[][]namesLocal,
    char[]startLoop,char[] loopBody,char[]endLoop,char[] uniq="")
{
    bool dataLoop=(pFlags&ParaFlags.DataLoop)!=0;
    auto pFlagsBase= (pFlags&ParaFlags.MaskBase);
    char[] visitorStr=`
    {
        `~iterName~` mainContext`~uniq~`;
        auto kEnd`~uniq~`=`~namesLocal[0]~`.kRange.kEnd;
        `~iterName~`* newK`~uniq~`;
        PIndex dummyP`~uniq~`;
        mainContext`~uniq~`.optimalBlockSize=optimalBlockSize`~uniq~`;
        mainContext`~uniq~`.context=&mainContext`~uniq~`;`;
    visitorStr~=startLoop;
    if (pFlagsBase!=ParaFlags.Sequential){
        visitorStr~=`
        Task("segArrayKLoop`~iterName~`",delegate void(){`;
    }
    visitorStr~=`
            for (KindIdx kIdx`~uniq~`=`~namesLocal[0]~`.kRange.kStart;kIdx`~uniq~`<kEnd`~uniq~`;++kIdx`~uniq~`){
                if (mainContext`~uniq~`.exception !is null) break;
                bool visitKind`~uniq~`=true;
                bool outOfRange`~uniq~`=false;
                if (newK`~uniq~` is null){
                    newK`~uniq~`=mainContext`~uniq~`.alloc();
                }
                newK`~uniq~`.kind=kIdx`~uniq~`;
                newK`~uniq~`.start=0;`;
    if (dataLoop) {
        visitorStr~=`
                auto aStruct=`~namesLocal[0]~`.arrayStruct;
                assert(aStruct.mKindDims[kIdx`~uniq~`-`~namesLocal[0]~`.kRange.kStart]!=0||(aStruct.kindStarts[kIdx`~uniq~`-aStruct.kRange.kStart+1]-aStruct.kindStarts[kIdx`~uniq~`-aStruct.kRange.kStart])==0);
                
                newK`~uniq~`.end=cast(ParticleIdx)
                    ((aStruct.kindStarts[kIdx`~uniq~`-aStruct.kRange.kStart+1]-aStruct.kindStarts[kIdx`~uniq~`-aStruct.kRange.kStart])/
                        max(1,aStruct.mKindDims[kIdx`~uniq~`-`~namesLocal[0]~`.kRange.kStart]));`;
    } else {
        visitorStr~=`
                newK`~uniq~`.end=
                    `~namesLocal[0]~`.arrayStruct.submapping.nLocalParticles(kIdx`~uniq~`);`;
    }
    visitorStr~=`
                auto submap`~uniq~`=`~namesLocal[0]~`.arrayStruct.submapping;
                if (submap`~uniq~`.kindStarts[kIdx`~uniq~`-submap`~uniq~`.lKRange.kStart]<
                    submap`~uniq~`.kindStarts[kIdx`~uniq~`-submap`~uniq~`.lKRange.kStart+1])
                {
                    newK`~uniq~`.pIndexPtrStart=submap`~uniq~`.ptrI(
                        LocalPIndex(newK`~uniq~`.kind,newK`~uniq~`.start));
                } else {
                    newK`~uniq~`.pIndexPtrStart=null;
                }`;
    foreach (i,name;namesLocal){
        visitorStr~=`
                if (kIdx`~uniq~` in `~name~`.kRange){
                    auto ik`~uniq~`=cast(size_t)(kIdx`~uniq~`-`~name~`.kRange.kStart);
                    if (`~name~`.mKindDims[ik`~uniq~`]!=0){
                        newK`~uniq~`.`~name~`PtrStart=cast(typeof(newK`~uniq~`.`~name~`PtrStart))(`~name~`.basePtr+`~name~`.kindOffsets[ik`~uniq~`]);
                        newK`~uniq~`.`~name~`MKindDims=`~name~`.mKindDims[ik`~uniq~`];
                        newK`~uniq~`.`~name~`ByteIncrements=`~name~`.kindByteIncrements[ik`~uniq~`];
                    } else {
                        outOfRange`~uniq~`=true;
                        newK`~uniq~`.`~name~`PtrStart=null;
                        newK`~uniq~`.`~name~`MKindDims=0;
                        newK`~uniq~`.`~name~`ByteIncrements=0;
                    }
                } else {
                    outOfRange`~uniq~`=true;
                    newK`~uniq~`.`~name~`PtrStart=null;
                    newK`~uniq~`.`~name~`MKindDims=0;
                    newK`~uniq~`.`~name~`ByteIncrements=0;
                }`;
    }
    visitorStr~="\n";
    visitorStr~=loopBody;
    visitorStr~="\n";
    visitorStr~=`
        }`;
    if (pFlagsBase!=ParaFlags.Sequential){
        visitorStr~=`
        }).autorelease.executeNow();`;
    }
    visitorStr~=`
        if (mainContext`~uniq~`.exception !is null) throw new Exception("exception in SegmentedArray loop",__FILE__,__LINE__,mainContext`~uniq~`.exception);
        auto freeL=mainContext`~uniq~`.next;
        while (freeL!is null){
            auto cNext=freeL.next;
            delete freeL;
            freeL=cNext;
        }
    `;
    visitorStr~=endLoop;
    visitorStr~=`
    }`;
    return visitorStr;
}

/// needs "optimalBlockSize"~uniq
/// public vars in kindVisit: struct iterName, "context"~uniq, "kIdx"~uniq, "outOfRange"~uniq, "visitKind"~uniq
/// public vars in pVisit: name~"Ptr"~uniq, "localPIndex"~uniq, "pIndexPtr"~uniq, "index"~uniq
/// there can be several versions of pVisit, if only one is given then that should be the generic version,
/// if two versions are given then the first is when there is one element per particle or per kind,
/// and the other is the generic, if 3 are given there is also one when they have the same size or 1
/// (i.e a single loop with increments 0 and 1 would work)
char[] segArrayMonoLoop(int pFlags,char[] iterName, char[][]namesLocal,
    char[] contextExtra,char[] uniq,char[]startLoop,char[]kindVisit,char[]endLoop,
    char[][]pVisitors)
{
    auto pFlagsBase=(pFlags &ParaFlags.MaskBase);
    char[] res="{\n";
    res~=segArrayContextStr(pFlags,iterName, namesLocal, contextExtra,uniq,
        pVisitors);
    char[] kVisit=kindVisit~`
            if (visitKind`~uniq~`){
                newK`~uniq~`.maxInnerDim=1;
                bool multiDim`~uniq~`=false;`;
    foreach (n;namesLocal){
                kVisit~=`
                if (newK`~uniq~`.`~n~`MKindDims!=newK`~uniq~`.maxInnerDim && newK`~uniq~`.`~n~`MKindDims!=0 && newK`~uniq~`.`~n~`MKindDims!=1){
                    if (newK`~uniq~`.maxInnerDim!=1){
                        multiDim`~uniq~`=true;
                    }
                    if (newK`~uniq~`.maxInnerDim<newK`~uniq~`.`~n~`MKindDims) newK`~uniq~`.maxInnerDim=newK`~uniq~`.`~n~`MKindDims;
                }`;
    }
    if (pVisitors.length==3) {
        if (pFlagsBase==ParaFlags.Sequential){
            kVisit~=`
                newK`~uniq~`.exec0();`;
        } else {
            kVisit~=`
                Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                    .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                newK`~uniq~`=null;`;
        }
    } else if (pVisitors.length==6){
        if (pFlagsBase==ParaFlags.Sequential){
            kVisit~=`
                if (newK`~uniq~`.maxInnerDim==1){
                    newK`~uniq~`.exec0();
                } else {
                    newK`~uniq~`.exec1();
                }
                newK`~uniq~`=null;`;
        } else {
            kVisit~=`
                if (newK`~uniq~`.maxInnerDim==1){
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                } else {
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec1)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                }
                newK`~uniq~`=null;`;
        }
    } else if (pVisitors.length==9){
        if (pFlagsBase==ParaFlags.Sequential){
            kVisit~=`
                if (!multiDim`~uniq~`){
                    if (newK`~uniq~`.maxInnerDim==1){
                        newK`~uniq~`.exec0();
                    } else {
                        newK`~uniq~`.exec1();
                    }
                } else {
                    newK`~uniq~`.exec2();
                }`;
        } else {
            kVisit~=`
                if (!multiDim`~uniq~`){
                    if (newK`~uniq~`.maxInnerDim==1){
                        Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                            .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                    } else {
                        Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec1)
                            .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                    }
                } else {
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec2)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                }`;
        }
    } else {
        assert(0,"invalid number of elements in pVisitors");
    }
    kVisit~=`
            }
    `;
    res~=segArrayKLoopStr(pFlags,iterName,namesLocal,startLoop,kVisit,endLoop,uniq);
    res~="\n}\n";
    return res;
}

/// needs optimalBlockSize_1, and optimalBlockSize_2 (for loop 1, loop 2)
/// loop 1 (external)
/// public vars in kindVisit1: struct iterName, "context_1", "kIdx_1", "outOfRange_1", "visitKind_1"
/// public vars in pVisit1: name~"Ptr_1", "localPIndex_1", "pIndexPtr_1", "index_1"
/// loop 1 (internal, the vaiables in loop1 and)
/// public vars in kindVisit2: struct iterName, "context_2", "kIdx_2", "outOfRange_2", "visitKind_2"
/// public vars in pVisit2: name~"Ptr_2", "localPIndex_2", "pIndexPtr_2", "index_2"
char[] segArrayBinLoop(int pFlags,char[] iterName, char[][]namesLocal1, char[] contextExtra1,
    char[]startLoop1,char[]kindVisit1,char[]endLoop1,
    char[][]pVisitors1,
    char[][]namesLocal2, char[] contextExtra2,char[]startLoop2,char[] kindVisit2,
    char[]endLoop2,char[][]pVisitors2)
{
    char[] preLocalLoop2=`
        ParticleIdx index_1=outerLoopContext.lIndex;
        LocalPIndex localPIndex_1=LocalPIndex(outerLoopContext.kind,index_1);
        PIndex *pIndexPtr_1=outerLoopContext.pIndexPtrStart+cast(size_t)index_1;`;
    foreach(name;namesLocal1){
        preLocalLoop2~=`
        auto `~name~`Ptr_1=cast(`~name~`.dtype *)(cast(size_t)outerLoopContext.`~name~`PtrStart
                +outerLoopContext.`~name~`ByteIncrements*cast(size_t)outerLoopContext.lIndex;
        `;
    }
    char[][] pVis2;
    for(int i=0;i<pVisitors2.length;++i){
        if (i%3==0){
            pVis2~=preLocalLoop2~pVisitors2[i];
        } else {
            pVis2~=pVisitors2[i];
        }
    }
    auto intC=segArrayContextStr(pFlags,iterName~"_2", namesLocal2, contextExtra2~"\n"~iterName~"_1 *outerLoopContext;\n","_2",
    pVis2);
    
    char[][] pVis1;
    for(int i=0;i<pVisitors1.length;++i){
        if (i%3==1){
            pVis1~=`
            auto innerContext=innerLoopStartContext.alloc();
            innerLoopStartContext.outerLoopContext=this;
            `~pVisitors1[i]~`
            Task("particleLoop`~iterName~`",&innerContext.exec0).appendOnFinish(&innerContext.giveBack).autorelease.submitYield();
            `;
        } else {
            pVis1~=pVisitors1[i];
        }
    }

    auto extC=segArrayContextStr(pFlags,iterName~"_1", namesLocal1,contextExtra1~"\n"~intC~"\n"~iterName~"_2 *innerLoopStartContext;\n","_1",
        pVis1);
    auto kLoopInt=segArrayKLoopStr(pFlags,iterName~"_1."~iterName~"_2", namesLocal2,startLoop2,kindVisit2~`
    if (visitKind_2){
        newK_1.innerLoopStartContext=newK_2;
        newK_2.outerLoopContext=newK_1;
        auto tmpK1=mainContext_1.alloc();
        (*tmpK1)=(*newK_1);
        Task("kindIntLoop`~iterName~`",&newK_1.exec0).appendOnFinish(&newK_1.giveBack).autorelease.submitYield();
        newK_1=tmpK1;
    }
    `,endLoop2,"_2");
    auto kLoopExt=segArrayKLoopStr(pFlags,iterName~"_1", namesLocal1,startLoop1,kindVisit1~`
    if (visitKind_1){
        Task("kindMainLoop`~iterName~`",delegate void(){
        `~kLoopInt~`
        }).autorelease.submitYield();
    }
    `,endLoop1,"_1");
    
    char[] res=`{
    struct `~iterName~"_2;\n";
    res~=extC;
    res~=kLoopExt;
    res~="\n}\n";
    return res;
}
