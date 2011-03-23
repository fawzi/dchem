/// particle indexes
module dchem.sys.PIndexes;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import blip.parallel.smp.WorkManager;
import blip.container.Cache;
import blip.container.Pool;
import blip.core.sync.Mutex;

version(LongIndexes){
    alias ulong idxType; /// type used for particle,... indexing
    typedef ulong ParticleIdx=0xFFFF_FFFF_FFFFUL; /// indexes a particle within a kind
} else {
    alias uint idxType; /// type used for particle,... indexing
    typedef uint ParticleIdx=0xFFFF_FFFF; /// indexes a particle within a kind
}

typedef ubyte LevelIdx=0xFF; /// indexes the level of the particle, level 255 is invalid
typedef ushort KindIdx=0xFFFF; /// indexes the kind of the particle, kind 65535 is invalid, kinds are global

static this(){
    {
        auto metaI=ClassMetaInfo.createForType!(ParticleIdx)("dchem.sys.PIndexes.ParticleIdx","index of a particle (within those of the same kind)");
        metaI.kind=TypeKind.CustomK;
        auto extH=new ExternalSerializationHandlers;
        metaI.externalHandlers=extH;
        extH.serialize=function void(Serializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ulong*)o);
        };
        extH.unserialize=function void(Unserializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ulong*)o);
        };
    }
    {
        auto metaI=ClassMetaInfo.createForType!(LevelIdx)("dchem.sys.PIndexes.LevelIdx","index of a level of the particles (0 are the atoms)");
        metaI.kind=TypeKind.CustomK;
        auto extH=new ExternalSerializationHandlers;
        metaI.externalHandlers=extH;
        extH.serialize=function void(Serializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ubyte*)o);
        };
        extH.unserialize=function void(Unserializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ubyte*)o);
        };
    }
    {
        auto metaI=ClassMetaInfo.createForType!(KindIdx)("dchem.sys.PIndexes.KindIdx","index identifying a kind of particles");
        metaI.kind=TypeKind.CustomK;
        auto extH=new ExternalSerializationHandlers;
        metaI.externalHandlers=extH;
        extH.serialize=function void(Serializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ushort*)o);
        };
        extH.unserialize=function void(Unserializer s,ClassMetaInfo mInfo,void* o){
            s.handlers.handle(*cast(ushort*)o);
        };
    }
}

/// a range of kinds
struct KindRange{
    KindIdx kStart; // start of the range (inclusive)
    KindIdx kEnd; // end of the range (exclusive)
    mixin(serializeSome("KindRange","A range of kinds (types of particles).","kStart|kEnd"));
    
    static KindRange opCall(KindIdx kStart, KindIdx kEnd){
        KindRange res;
        res.kStart=kStart;
        res.kEnd=kEnd;
        if (res.kEnd<res.kStart){
            res.kEnd=res.kStart; // avoid?
        }
        return res;
    }
    /// a pseudo range that represents an infinite range (not accepted everywhere)
    static KindRange all(){
        KindRange res;
        res.kStart=cast(KindIdx)0;
        return res;
    }
    /// returns the intersection of this range with kr
    KindRange intersect(KindRange kr){
        KindRange res;
        res.kStart=((kStart>kr.kStart)?kStart:kr.kStart);
        res.kEnd=((kEnd<kr.kEnd)?kEnd:kr.kEnd);
        if (res.kEnd<res.kStart){
            res.kEnd=res.kStart;
        }
        return res;
    }
    /// returns the convex hull of this range and kr (empy ranges are ignored)
    KindRange convexHull(KindRange kr){
        if (kr.kStart>=kr.kEnd){
            return *this;
        } else if (kStart>=kEnd){
            return kr;
        } else {
            KindRange res;
            res.kStart=((kStart<kr.kStart)?kStart:kr.kStart);
            res.kEnd=((kEnd>kr.kEnd)?kEnd:kr.kEnd);
            return res;
        }
    }
    bool opIn_r(KindIdx i){
        return i>=kStart && i<kEnd;
    }
    bool opIn(KindRange kr){
        return kr.kStart<=kStart && kr.kEnd>=kEnd;
    }
    int opApply(int delegate(ref KindIdx)loopOp){
        for (KindIdx k=kStart;k<kEnd;++k){
            auto res=loopOp(k);
            if (res!=0) return res;
        }
        return 0;
    }
    /// parallel loop on this kind range
    struct PLoop{
        int delegate(ref KindIdx) loopOp;
        Exception exception;
        int optimalBlockSize=1;
        int res=0;
        KindRange kr;
        
        static PLoop opCall(KindRange k,int blockSize=1){
            PLoop res;
            assert(blockSize>0,"blockSize cannot be 0");
            res.kr=k;
            res.optimalBlockSize=blockSize;
            return res;
        }

        struct LoopK{
            PLoop *ctx;
            KindRange kr;
            PoolI!(LoopK*) pool;
            LoopK *next;
            
            static Mutex gLock;
            static CachedPool!(LoopK*) gPool;
            static size_t nPool;
            static this(){
                gLock=new Mutex();
            }
            static void addGPool(){
                synchronized(gLock){
                    if (nPool==0){
                        assert(gPool is null);
                        gPool=cachedPoolNext(delegate LoopK*(PoolI!(LoopK*)p){
                            auto res=new LoopK;
                            res.pool=p;
                            return res;
                        });
                    }
                    nPool+=1;
                }
            }
            static void rmGPool(){
                synchronized(gLock){
                    if (nPool==0){
                        throw new Exception("mismatched rmGPool",__FILE__,__LINE__);
                    }
                    --nPool;
                    if (nPool==0){
                        assert(gPool!is null);
                        gPool.stopCaching();
                        gPool=null;
                    }
                }
            }
            static void gSync(){
                gLock.lock(); gLock.unlock();
            }
            void doOp(){
                auto bSize=ctx.optimalBlockSize;
                if (ctx.exception!is null || ctx.res!=0) return;
                if (kr.kEnd>kr.kStart+bSize*3/2){
                    if (gPool is null) gSync();
                    assert(gPool!is null,"gPool is null (missin addGPool?)");
                    auto k1=gPool.getObj();
                    auto k2=gPool.getObj();
                    int mid=(kr.kEnd-kr.kStart)/2;
                    if (mid>bSize){
                        mid=((mid+bSize-1)/bSize)*bSize;
                    }
                    mid+=cast(int)kr.kStart;
                    k1.kr.kStart=kr.kStart;
                    k1.kr.kEnd=cast(KindIdx)mid;
                    k2.kr.kStart=cast(KindIdx)mid;
                    k2.kr.kEnd=kr.kEnd;
                    Task("KindRangePLoopSub1",&k1.doOp).appendOnFinish(&k1.giveBack).autorelease.submit();
                    Task("KindRangePLoopSub2",&k2.doOp).appendOnFinish(&k2.giveBack).autorelease.submit();
                } else {
                    try{
                        auto loopB=ctx.loopOp;
                        for (KindIdx kIdx=kr.kStart;kIdx<kr.kEnd;++kIdx){
                            auto resAtt=loopB(kIdx);
                            if (resAtt!=0){
                                ctx.res=resAtt;
                                return;
                            }
                        }
                    } catch (Exception e){
                        ctx.exception=e;
                    }
                }
            }
            
            void giveBack(){
                if (pool!is null){
                    pool.giveBack(this);
                } else {
                    delete this;
                }
            }
        }
        
        int opApply(int delegate(ref KindIdx)loopOp){
            version(PLoopKinds){
                this.loopOp=loopOp;
                if (kr.kStart>=kr.kEnd) return 0;
                LoopK.addGPool();
                auto mainLooper=LoopK.gPool.getObj();
                mainLooper.ctx=this;
                mainLooper.kr=kr;
                Task("KindRangePLoop",&mainLooper.doOp).appendOnFinish(&mainLooper.giveBack).autorelease.executeNow();
                LoopK.rmGPool();
                if (exception!is null) throw new Exception("Exception in parallel kind loop",__FILE__,__LINE__,exception);
                return res;
            } else {
                for (auto ik=kr.kStart;ik<kr.kEnd;++ik){
                    auto res=loopOp(ik);
                    if (res!=0) return res;
                }
                return 0;
            }
        }
    }
    /// parallel loop (all kinds in parallel)
    PLoop pLoop(){
        return PLoop(*this);
    }
    /// true if the current KindRange is a dummy
    bool isDummy(){
        return kStart==KindIdx.init && kEnd==KindIdx.init;
    }
    mixin printOut!();
    /// true if the current KindRange is invalid or dummy
    bool valid(){
        return kStart<=kEnd && kEnd<KindIdx.init;
    }
    size_t length(){
        assert(kEnd!=KindIdx.init);
        return kEnd-kStart;
    }
}

/// particle index, stores particle level,particle kind and particle position in one single long
/// at the moment there is one unused byte that is always 0, it could be used to increase the range
/// of the particles (for example)
struct PIndex{
    ulong data=0xFFFF_FFFF_FFFF_FFFFUL;
    static if(ParticleIdx.sizeof==4){
        enum:ulong{ ParticleMask=0x0000_0000_FFFF_FFFFUL }
    } else {
        enum:ulong{ ParticleMask=0x0000_FFFF_FFFF_FFFFUL }
    }
    
    static PIndex opCall()(){
        PIndex res;
        return res;
    }
    static PIndex opCall(T,U)(T k, U p){
        PIndex res;
        static if (!is(T==KindIdx)){
            static assert(!(is(T==ParticleIdx)||is(T==LevelIdx)),"invalid type "~T.stringof);
            assert((cast(ulong)k)<=0xFFFF,"kind out of range");
        }
        static if (!is(U==ParticleIdx)){
            static assert(!(is(U==KindIdx)||is(U==LevelIdx)),"invalid type "~U.stringof);
            assert((cast(ulong)p)<=cast(ulong)ParticleIdx.init,"particle out of range");
        }
        res.data=((cast(ulong)k)<<48)|(cast(ulong)p);
        return res;
    }
    static PIndex opCall(T)(T p){
        PIndex res;
        static if (is(T==ulong)){
            res.data=p;
        } else static if (is(T==PIndex)){
            return p;
        } else static if (is(T==LocalPIndex)){
            res.data=p.data;
        } else {
            static assert(0,"unsupported type "~T.stringof~" in PIndex.opCall");
        }
        return res;
    }
    /// returns the kind of the particle
    /// kinds are supposed to be global, and ordered in level sequence without gaps
    KindIdx kind(){
        return cast(KindIdx)((data>>48)&0xFFFFUL);
    }
    /// sets the kind of the particle
    void kindT(T)(T newK){
        static if (!is(T==KindIdx)){
            static assert(!(is(T==ParticleIdx)||is(T==LevelIdx)),"invalid type "~T.stringof);
            assert((cast(ulong)newK)<=0xFFFF,"kind out of range");
        }
        data=(data & 0x0000_FFFF_FFFF_FFFFUL)|((cast(ulong)newK)<<48);
    }
    alias kindT!(KindIdx) kind;
    alias kindT!(short) kind;
    alias kindT!(ushort) kind;
    alias kindT!(int) kind;
    alias kindT!(uint) kind;
    alias kindT!(long) kind;
    alias kindT!(ulong) kind;
    ParticleIdx particle(){
        return cast(ParticleIdx)(data & 0xFFFF_FFFF_FFFFUL);
    }
    void particleT(T)(T pNew){
        static if (!is(T==ParticleIdx)){
            static assert(!(is(T==KindIdx)||is(T==LevelIdx)),"invalid type "~T.stringof);
            assert((cast(ulong)pNew)<=cast(ulong)ParticleIdx.init,"particle out of range");
        }
        data=(data & 0xFFFF_0000_0000_0000UL)|(cast(ulong) pNew);
    }
    alias particleT!(ParticleIdx) particle;
    alias particleT!(int) particle;
    alias particleT!(uint) particle;
    alias particleT!(long) particle;
    alias particleT!(ulong) particle;
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.sys.PIndexes.PIndex",
            "a global index of a particle");
        metaI.addFieldOfType!(ushort)("kind","the kind of the particle");
        metaI.addFieldOfType!(uint)("particle","the index of the particle");
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    void serialize(Serializer s){
        ushort us=kind;
        s.field(metaI[0],us);
        uint ui=particle;
        s.field(metaI[1],ui);
    }
    void unserialize(Unserializer s){
        ushort us=kind;
        s.field(metaI[0],us);
        kind=cast(KindIdx)us;
        uint ui=particle;
        s.field(metaI[1],ui);
        particle=cast(ParticleIdx)ui;
    }
    /// quick check to see if the particle is a dummy particle
    bool dummy(){
        return data==0xFFFF_FFFF_FFFF_FFFFUL;
    }
    /// slower check that looks if any component is invalid
    bool valid(){
        return ((data & ParticleMask) != ParticleMask) &&
            ((data & 0xFFFF_0000_0000_0000UL) != 0xFFFF_0000_0000_0000UL);
    }
    // default implementations are ok
    // int opCmp(PIndex o)
    // equals_t    opEquals(Object o);
    
    /// increments the particle part
    PIndex opAddAssign(int i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// ditto
    PIndex opAddAssign(ulong i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// ditto
    PIndex opAddAssign(ParticleIdx i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// increments the kind part
    PIndex opAddAssign(KindIdx i){
        assert((data & 0xFFFF_0000_0000_0000UL) != 0xFFFF_0000_0000_0000UL,"kind increment invalid particle kind");
        data+=(cast(ulong)i)<<48;
        return *this;
    }
    /// adds a value (only on the right side, add also the left??)
    PIndex opAdd(T)(T v){
        PIndex res=*this;
        return (res+=v);
    }
    /// compares two PIndex
    int opCmp(PIndex p2){
        return ((this.data<p2.data)?-1:((this.data==p2.data)?0:1));
    }
    mixin printOut!();
}

/// an index local to a subset of atoms (clearly not transferrable between different subsets)
// typedef is too limited
struct LocalPIndex{
    ulong data=0xFFFF_FFFF_FFFF_FFFFUL;
    static if(ParticleIdx.sizeof==4){
        enum:ulong{ ParticleMask=0x0000_0000_FFFF_FFFFUL }
    } else {
        enum:ulong{ ParticleMask=0x0000_FFFF_FFFF_FFFFUL }
    }
    
    static LocalPIndex opCall()(){
        LocalPIndex res;
        return res;
    }
    static LocalPIndex opCall(T,U)(T k, U p){
        LocalPIndex res;
        static if (!is(T==KindIdx)){
            static assert(!(is(T==ParticleIdx)||is(T==LevelIdx)),"invalid type "~T.stringof);
            assert((cast(ulong)k)<=0xFFFF,"kind out of range");
        }
        static if (!is(U==ParticleIdx)){
            static assert(!(is(U==KindIdx)||is(U==LevelIdx)),"invalid type "~U.stringof);
            assert((cast(ulong)p)<=cast(ulong)ParticleIdx.init,"particle out of range");
        }
        res.data=((cast(ulong)k)<<48)|(cast(ulong)p);
        return res;
    }
    static LocalPIndex opCall(T)(T p){
        LocalPIndex res;
        static if (is(T==ulong)){
            res.data=p;
        } else static if (is(T==LocalPIndex)){
            return p;
        } else static if (is(T==PIndex)){
            res.data=p.data;
        } else {
            static assert(0,"unsupported type "~T.stringof~" in PIndex.opCall");
        }
        return res;
    }
    /// returns the kind of the particle
    /// kinds are supposed to be global, and ordered in level sequence without gaps
    KindIdx kind(){
        return cast(KindIdx)((data>>48)&0xFFFFUL);
    }
    /// sets the kind of the particle
    void kindT(T)(T newK){
        static if (!is(T==KindIdx)){
            static assert(!(is(T==ParticleIdx)||is(T==LevelIdx)),"invalid type"~T.stringof);
            assert((cast(ulong)newK)<=0xFFFF,"kind out of range"); // remove check?
        }
        data=(data & 0x0000_FFFF_FFFF_FFFFUL)|((cast(ulong)newK)<<48);
    }
    alias kindT!(KindIdx) kind;
    alias kindT!(short) kind;
    alias kindT!(ushort) kind;
    alias kindT!(int) kind;
    alias kindT!(uint) kind;
    alias kindT!(long) kind;
    alias kindT!(ulong) kind;
    ParticleIdx particle(){
        return cast(ParticleIdx)(data & ParticleMask);
    }
    void particleT(T)(T pNew){
        static if (!is(T==ParticleIdx)){
            static assert(!(is(T==KindIdx)||is(T==LevelIdx)),"invalid type "~T.stringof);
            assert((cast(ulong)pNew)<=cast(ulong)ParticleIdx.init,"particle out of range"); // remove check?
        }
        data=(data & 0xFFFF_0000_0000_0000UL)|(cast(ulong) pNew);
    }
    alias particleT!(ParticleIdx) particle;
    alias particleT!(int) particle;
    alias particleT!(uint) particle;
    alias particleT!(long) particle;
    alias particleT!(ulong) particle;
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.sys.PIndexes.LocalPIndex",
            "a local index of a particle");
        metaI.addFieldOfType!(ushort)("kind","the kind of the particle");
        metaI.addFieldOfType!(uint)("particle","the index of the particle");
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    void serialize(Serializer s){
        ushort us=kind;
        s.field(metaI[0],us);
        uint ui=particle;
        s.field(metaI[1],ui);
    }
    void unserialize(Unserializer s){
        ushort us=kind;
        s.field(metaI[0],us);
        kind=cast(KindIdx)us;
        uint ui=particle;
        s.field(metaI[1],ui);
        particle=cast(ParticleIdx)ui;
    }
    /// quick check to see if the particle is a dummy particle
    bool dummy(){
        return data==0xFFFF_FFFF_FFFF_FFFFUL;
    }
    /// slower check that looks if any component is invalid
    bool valid(){
        return ((data & ParticleMask) != ParticleMask) &&
            ((data & 0xFFFF_0000_0000_0000UL) != 0xFFFF_0000_0000_0000UL);
    }
    // default implementations are ok
    // int opCmp(LocalPIndex o)
    // equals_t    opEquals(Object o);
    
    /// increments the particle part
    LocalPIndex opAddAssign(int i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// ditto
    LocalPIndex opAddAssign(ulong i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// ditto
    LocalPIndex opAddAssign(ParticleIdx i){
        assert((data & ParticleMask) != ParticleMask,"increment invalid particle");
        data+=i;
        return *this;
    }
    /// increments the kind part
    LocalPIndex opAddAssign(KindIdx i){
        assert((data & 0xFFFF_0000_0000_0000UL) != 0xFFFF_0000_0000_0000UL,"kind increment invalid particle kind");
        data+=(cast(ulong)i)<<48;
        return *this;
    }
    /// adds a value (only on the right side, add also the left??)
    LocalPIndex opAdd(T)(T v){
        LocalPIndex res=*this;
        return (res+=v);
    }
    /// compares two LocalPIndex
    int opCmp(LocalPIndex p2){
        return ((this.data<p2.data)?-1:((this.data==p2.data)?0:1));
    }
    mixin printOut!();
}

