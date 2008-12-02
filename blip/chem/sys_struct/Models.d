module blip.chem.sysStruct.Models;

alias int idxType; /// type used for particle,... indexing

/// forward iterator interface
interface FIterator(T){
    T next();
    bool atEnd();
    int opApply(int delegate(T x) dlg); // might be parallel
    int opApply(int delegate(idxType i,T x) dlg); // might be parallel
}

/// Particle
/// minimal information on a particle, this is a structure, and should be kept small 
struct Particle
{
    idxType globalIdx; /// index of the particle (its position in the global array)
    idxType kindIdx; /// index of the kind of the current particle
    idxType inKindIdx; /// index of the particle within the kind
    idxType superPartIdx; /// index of the super particle of the current particle (-1 if none)
    idxType inSuperPartIdx; /// index of the current particle in its superParticle
    Particle *superParticle(SysStruct s);
}

/// represents the mappings at a given level
interface LevelMappings{
    int level();
    FIterator!(idxType) subParticleIdxs(idxType pKind, idxType localParticleIdx);
    FIterator!(idxType) bottomParticlesIdxs(idxType pKind, idxType localParticleIdx);
    FIterator!(Particle) subParticles(idxType pKind, idxType localParticleIdx);
    FIterator!(Particle) bottomParticles(idxType pKind, idxType localParticleIdx);
}

/// represent a group of particles with the same kind
interface ParticleKindKind{
    char[] particleKindName();
    int level();
    bool fixedLayout();
    PropertyKind[] propertyKinds();
    PropertyKind getPropertyKind(char[] propertyName)
    ParticleKind[] subParticleKinds;
    DynPPropKind[] dynPPropKinds;
    Constraints[] constaints;
    void serialize(Serializer s);
    ParticleKind dup();
}


/// represent a group of particles with the same kind
interface ParticleKind{
    char[] particleKindName();
    int level();
    int pKindIdx();
    GenProperty[] properties();
    GenProperty getProperty(char[] propertyName)
    FIterator!(idxType) particles();
    FIterator!(idxType) subParticleIdxs(idxType localParticleIdx);
    FIterator!(idxType) bottomParticlesIdxs(idxType localParticleIdx);
    FIterator!(Particle) subParticles(idxType localParticleIdx);
    FIterator!(Particle) bottomParticles(idxType localParticleIdx);
    void serialize(Serializer s);
    ParticleKind dup();
    ParticleKindKind kind();
}

/// Meta informations, properties of a property
interface PropertyKind{
    enum Distribution:int{
        Replicated,
        ParticleLocal
    }
    enum DetailLevel:int{
        KindLevel,
        ParticleLevel
    }
    enum Storage:int{
        NArrayReal0DT,
        NArrayReal1DT,
        NArrayReal2DT,
        NArrayReal3DT,
        VariantT
    }
    Distribution distribution();
    DetailLevel detailLevel();
    Storage storageLevel();
    char[] propertyName();
    int propertyIdx();
    ParticleKind particleKind();
}

/// generic property interface
interface GenProperty {
    PropertyKind kind(); /// kind of this property
    Variant storage(); /// the data
    void serialize(Serializer s); /// serializes the data
    GenProperty dup(); /// duplicates the property
}

/// properties local to a particle
interface ParticleProperty(T): GenProperty{
    T opIndex(idxType i);
    void opIndexAssign(idxType i,T val);
    void value(T val);
}

/// property local to a kind
interface PKindProperty(T): GenProperty{
    T value();
    void value(T val);
    T opIndex(idxType i);
}

