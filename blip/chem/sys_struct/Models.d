module blip.chem.sysStruct.Models;

alias int idxType; /// type used for particle,... indexing

/// forward iterator interface
interface FIteratorI(T){
    T next();
    bool atEnd();
    int opApply(int delegate(T x) dlg);
    int opApply(int delegate(idxType i,T x) dlg);
    FIteratorI(T) parallelIterator();
}

/// basic interface for objects
interface DuplicableI{
    Duplicable dup();
}

/// basic interface for objects
interface DeepDuplicableI{
    Duplicable deepdup();
}

interface CopiableOjectI : BasicObjectI,DuplicableI,DeepDuplicableI { }

/// represent a kind of particle (has no direct pointer to actual particles)
interface ParticleKindKind:CopiableOjectI{
    char[] particleKindName();
    int level();
    bool fixedLayout();
    PropertyKind[] propertyKinds();
    PropertyKind getPropertyKind(char[] propertyName)
    ConstraintsKinds[] constaints;
    void serialize(Serializer s);
}


/// represent a group of particles with the same kind
interface ParticleKind: CopiableOjectI{
    char[] particleKindName();
    int level();
    int pKindIdx();
    GenProperty[] properties();
    GenProperty propertyNamed(char[] propertyName);
    FIteratorI!(idxType) particles();
    FIteratorI!(idxType) subParticleIdxs(idxType localParticleIdx);
    FIteratorI!(idxType) bottomParticlesIdxs(idxType localParticleIdx);
    FIteratorI!(Particle) subParticles(idxType localParticleIdx);
    FIteratorI!(Particle) bottomParticles(idxType localParticleIdx);
    void serialize(Serializer s);
    ParticleKindKind kind();
    bool generateLevelKindIndexes(); /// if the level info should keep indexes about the particles with this kind
    bool keepSubparticleIndexes(); /// if the level info should keep indexes about subparticles
}

/// Meta informations, properties of a property
interface PropertyKind: CopiableOjectI{
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
interface GenProperty:CopiableOjectI {
    PropertyKind kind(); /// kind of this property
    Variant storage(); /// the data
    void serialize(Serializer s); /// serializes the data
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


/// Particle
/// minimal information on a particle, this is a structure, and should be kept small 
struct Particle
{
    idxType globalIdx; /// index of the particle (its position in the global array)
    int level; /// level of the current particle
    idxType kindIdx; /// index of the kind of the current particle
    idxType inKindIdx; /// index of the particle within the kind
    idxType superParticleIdx; /// index of the super particle of the current particle (-1 if none)
    idxType inSuperParticleIdx; /// index of the current particle in its superParticle
    Particle *superParticle(SysStruct s); /// return a pointer to the superParticle of this particle
    ParticleKind kind(SysStruct s); /// returns the kind of the current particle
}

/// represents the mappings at a given level
class LevelMappings: CopiableOjectI {
    this(){}
    NArray!(idxType,1)[] kinds2particles;
    NArray!(idxType,1)[] superParticles2particles;
    void updateWithParticles(ParticleSys psys,int level){
        kinds2particles.length=kinds.length;
        superParticles2particles=superParticles.length;
        idxType[] nPartPerKind=new idxType[kinds.length];
        idxType[] nPartPerSuperP=new idxType[superParticles.length];
        for(idxType i=0;i!=particles.length;++i){
            particles[i]=i;
            nPartPerKind[particles[i].kindIdx]+=1;
            nPartPerSuperP[particles[i].superParticleIdx]+=1;
        }
        foreach(i,k;kinds){
            if (!k.keepSubparticleIndexes()) nPartPerKind[i]=0;
        }
        foreach(i,sp;superParticles){
            if (!sp.keepSubparticleIndexes()) nPartPerSuperP[i]=0;
        }
        for(idxType i=0;i!=particles.length;++i){
            particles[i]=i;
            nPartPerKind[particles[i].kindIdx]+=1;
            nPartPerSuperP[particles[i].superParticleIdx]+=1;
        }
        
    }
    FIteratorI!(idxType) subParticleIdxs(idxType pKind, idxType localParticleIdx);
    FIteratorI!(idxType) bottomParticlesIdxs(idxType pKind, idxType localParticleIdx);
    FIteratorI!(Particle) subParticles(idxType pKind, idxType localParticleIdx);
    FIteratorI!(Particle) bottomParticles(idxType pKind, idxType localParticleIdx);
}

/// represent a system of particles
class ParticleSys
{
    ParticleKinds[][] pKinds;
    Particle[][] particles; // the particles are mallocated, and freed by this class
    LevelMappings[] pLevels;
    DynSysPropKind[] sysProperties();
    ParticleSys dup();
    ParticleSys deepDup();
    FIteratorI!(Particle *) particleIterator(int level);
}

// BulkArray for particles: 1D, malloc, slices, indexing/itarations return pointers, parallel iterators