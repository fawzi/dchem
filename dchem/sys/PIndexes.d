/// particle indexes
module dchem.sys.PIndexes;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import blip.parallel.smp.WorkManager;

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
        auto metaI=ClassMetaInfo.createForType!(ParticleIdx)("dchem.sys.PIndexes.ParticleIdx");
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
        auto metaI=ClassMetaInfo.createForType!(LevelIdx)("dchem.sys.PIndexes.LevelIdx");
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
        auto metaI=ClassMetaInfo.createForType!(KindIdx)("dchem.sys.PIndexes.KindIdx");
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
    mixin(serializeSome("dchem.sys.KindRange","kStart|kEnd"));
    
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
    int opApply(int delegate(KindIdx)loopOp){
        for (KindIdx k=kStart;k<kEnd;++k){
            auto res=loopOp(k);
            if (res!=0) return res;
        }
        return 0;
    }
    /// parallel loop on this kind range
    struct PLoop{
        KindRange kr;
        static PLoop opCall(KindRange k){
            PLoop res;
            res.kr=k;
            return res;
        }
        int opApply(int delegate(KindIdx)loopOp){
            if (kr.kStart>=kr.kEnd) return 0;
            Exception exception;
            int res=0;
            struct LoopK{
                int delegate(KindIdx) inOp;
                KindIdx k;
                Exception* ep;
                int *resPtr;
                void doOp(){
                    if ((*ep) is null && *resPtr==0){
                        try{
                            auto resAtt=inOp(k);
                            if (resAtt!=0) *resPtr=resAtt;
                        } catch (Exception e){
                            *ep=e;
                        }
                    }
                }
                void seppuku(){
                    delete this;
                }
            }
            Task("KindRangePLoop",delegate void(){
                for (KindIdx k=kr.kStart;k<kr.kEnd;++k){
                    auto lNow=new LoopK;
                    lNow.inOp=loopOp;
                    lNow.k=k;
                    lNow.ep=&exception;
                    lNow.resPtr=&res;
                    Task("KindRangePLoopIter",&lNow.doOp).appendOnFinish(&lNow.seppuku).autorelease.submitYield();
                }
            }).autorelease.executeNow();
            if (exception!is null) throw new Exception("Exception in parallel kind loop",__FILE__,__LINE__,exception);
            return res;
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
        metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.sys.PIndexes.PIndex");
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
        metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.sys.PIndexes.LocalPIndex");
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

