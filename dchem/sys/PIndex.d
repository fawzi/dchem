/// particle indexes
module dchem.sys.ParticleIndexes;

version(LongIndexes){
    alias ulong idxType; /// type used for particle,... indexing
    typedef ulong ParticleIdx=0xFFFF_FFFF_FFFFUL; /// indexes a particle within a kind
} else {
    alias uint idxType; /// type used for particle,... indexing
    typedef uint ParticleIdx=0xFFFF_FFFF; /// indexes a particle within a kind
}

typedef ubyte LevelIdx=0xFF; /// indexes the level of the particle, level 255 is invalid
typedef ushort KindIdx=0xFFFF; /// indexes the kind of the particle, kind 65535 is invalid, kinds are global

/// a range of kinds
struct KindRange{
    KindIdx kStart; // start of the range (inclusive)
    KindIdx kEnd; // end of the range (exclusive)
    static KindRange opCall(KindIdx kStart, KindIdx kEnd){
        KindRange res;
        res.kStart=kStart;
        res.kEnd=kEnd;
        if (res.kEnd<res.kStart){
            res.kEnd=res.kStart; // avoid?
        }
        return res;
    }
    KindRange intersect(KindRange kr){
        KindRange res;
        res.kStart=((kStart>kr.kStart)?kStart:kEnd);
        res.kEnd=((kEnd<kr.kEnd)?kEnd:kr.kEnd);
        if (res.kEnd<res.kStart){
            res.kEnd=res.kStart; // avoid?
        }
        return res;
    }
    bool opIn(KindIdx i){
        return i>=kStart && i<kEnd;
    }
    bool opIn(KindRange kr){
        return kr.kStart<=kStart && kr.kEnd>=kEnd;
    }
    int opApply(int delegate(KindIdx)loopOp){
        for (KindIdx k=kStart;k<kEnd;++k){
            auto res=loopOp(k);
            ïf (res!=0) break;
        }
    }
    bool dummy(){
        return kStart==KindIdx() && kEnd==KindIdx();
    }
}

/// particle index, stores particle level,particle kind and particle position in one single long
/// at the moment there is one unused byte that is always 0, it could be used to increase the range
/// of the particles (for example)
struct PIndex{
    ulong data=0xFFFF_FF00_FFFF_FFFFUL;
    
    static PIndex opCall(){
        PIndex res;
        return res;
    }
    static PIndex opCall(KindIdx k, ParticleIdx p){
        PIndex res;
        res.data=((cast(ulong)k)<<48)|(cast(ulong)p);
        return res;
    }
    /// returns the kind of the particle
    /// kinds are supposed to be global, and ordered in level sequence without gaps
    KindIdx kind(){
        return cast(KindIdx)((data>>48)&0xFFFFUL);
    }
    /// sets the kind of the particle
    void kind(KindIdx newK){
        data=(data & 0x0000_FFFF_FFFF_FFFFUL)|((cast(ulong)newK)<<40);
    }
    ParticleIdx particle(){
        return cast(ParticleIdx)(data & 0xFFFF_FFFF_FFFFUL);
    }
    void particle(ParticleIdx pNew){
        data=(data & 0xFFFF_0000_0000_0000UL)|(cast(ulong) pNew);
    }
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(this))("PIndex");
        metaI.addFieldOfType!(ushort)("kind","the kind of the particle");
        metaI.addFieldOfType!(uint)("particle","the index of the particle");
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    void serial(S)(S s){
        s.field(metaI[0],cast(ushort)kind);
        s.field(metaI[1],cast(uint)particle);
    }
    void serialize(Serializer s){
        serial(s);
    }
    void unserialize(Unserializer s){
        serial(s);
    }
    /// quick check to see if the particle is a dummy particle
    bool dummy(){
        return data==0xFFFF_FFFF_FFFF_FFFFUL;
    }
    /// slower check that looks if any component is invalid
    bool valid(){
        return ((data & 0x0000_FFFF_FFFF_FFFFUL) != 0x0000_FFFF_FFFF_FFFFUL) &&
            ((data & 0xFFFF_0000_0000_0000UL) != 0xFFFF_0000_0000_0000UL);
    }
    // default implementations are ok
    // int opCmp(PIndex o)
    // equals_t    opEquals(Object o);
}

typedef PIndex LocalPIndex; /// an index local to a subset of atoms (clearly not transferrable between different subsets)
