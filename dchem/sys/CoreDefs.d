module dchem.sys.CoreDefs;
import dchem.sys.Requests;
import dchem.PIndexes;

/// various levels of duplication
enum PSCopyDepthLevel{
    None, /// does not copy
    PSysLevel, /// shallowest level, just duplicates the ParticleSys
    DynProperties, /// copies also normal particle properties (position,...)
    MostProperties, /// most properties
    SysStruct, /// copies particles and particle kinds
    DeepProperties, /// all properties
    All /// copies all
}

/// represents a constraints
interface Constraint{
    void applyR(); // args?
    void applyDR();
    void applyDDR();
    FIteratorI!(PIndex)particlesInvolved();
}

/// represent a group of particles with the same kind
class ParticleKind: CopiableOjectI,Serializable{
    char[] _name;
    char[] name(){ return _name; }
    void name(char[] nName){ _name=nName; }
    LevelIdx _level;
    LevelIdx level(){ return _level; }
    void level(LevelIdx l){ _level=l; }
    KindIdx pKind;
    size_t position;
    size_t degreesOfFreedom;
    size_t orientation;
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(this))("ParticleKind");
        metaI.addFieldOfType!(char[])("name","name of this particle kind",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(ubyte)("level","level of the particle",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(ushort)("pKind","kind index");
        metaI.addFieldOfType!(size_t)("position","number of 3D space elements",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(size_t)("degreesOfFreedom","number of 3D space elements",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(size_t)("orientation","number of 3D space elements",
            SerializationLevel.debugLevel);
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    
    void serial(S s){
        s.field(metaI[0],_name);
        s.field(metaI[1],_level);
        s.field(metaI[2],_pKind);
        s.field(metaI[3],position);
        s.field(metaI[1],degreesOfFreedom);
        s.field(metaI[1],orientation);
    }
    
    void preSerialize(Serializer s){}
    void postSerialize(Serializer s){}
    Serializable preUnserialize(Unserializer s){ return this; }
    Serializable postUnserialize(Unserializer s){ return this; }
    
    void serialize(Serializer s){
        serial(s);
    }
    void unserialize(Unserializer s){
        serial(s);
    }
    
    // callbacks
    /// system structure changed (particle added/removed, kinds added/removed)
    void sysStructChanged(ParticleSys p){}
    /// position of particles changed
    void positionsChanged(ParticleSys p){}
    /// cell changed
    void cellChanged(ParticleSys p){}
    /// properties to be calculated did change
    void propertiesRequestedChanged(ParticleSys p){}
    /// properties were allocated
    void propertiesAllocated(ParticleSys p){}
    /// will calculate the properties
    void willCalculate(ParticleSys p){}
    /// did calculate the properties
    void didCalculate(ParticleSys p){}
    /// did read the properties (no need to keep them in memory)
    void didReadProperties{}
    /// request to try to reduce memory usage
    void minimizeMemory(ParticleSys p){}
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars{
    // cell
    Cell cell;

    // position in 3D space (used for neigh lists/hierarchical partitioning, screening)
    SegmentedArray!(float[3]) spos;

    // position in 3D space
    SegmentedArray!(Real[3]) pos;
    SegmentedArray!(Real[3]) dpos;
    SegmentedArray!(Real[3]) ddpos;

    // orientation (quaternions)
    SegmentedArray!(Real[4]) orient;
    SegmentedArray!(Real[4]) dorient;
    SegmentedArray!(Real[4]) ddorient;

    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(Real) dof;
    SegmentedArray!(Real) ddof;
    SegmentedArray!(Real) dddof;
    
    // constraints
    SegmentedArray!(Constraints) constraints;
    
    DynamicsVars dup(PSCopyDepthLevel level){
        if (level>=PSCopyDepthLevel.DynamicsProperties){
            DynamicsVars res;
            res.cell=cell.dup;
            res.spos=spos.dup;
            res.pos=pos.dup;
            res.dpos=dpos.dup;
            res.ddpos=ddpos.dup;
            res.orient=orient.dup;
            res.dorient=orient.dup;
            res.ddorient=orient.dup;
            res.dof=dof.dup;
            res.ddof=ddof.dup;
            res.dddof=dddof.dup;
            return res;
        }
        return *this;
    }
    
    mixin(serializeSome("dchem.sys.DynamicsVars","cell|spos|pos|dpos|ddpos|orient|dorient|"))
}

/// represent the structure of a system of particles
class SysStruct, CopiableOjectI,Serializable
{
    SubMapping fullSystem;
    SegmentedArray!(PIndex) particles;
    SegmentedArray!(PIndex) superParticle;
    SegmentedArray!(size_t) subParticleIdxs;
    SegmentedArray!(PKinds) particleKinds;
    
    this(){ }
    this(SubMapping fullSystem,SegmentedArray!(PIndex) particles,
        SegmentedArray!(PIndex) superParticle,SegmentedArray!(size_t) subParticleIdxs,
        SegmentedArray!(PKinds) particleKinds)
    {
        this.fullSystem=fullSystem;
        this.particles=particles;
        this.superParticle=superParticle;
        this.subParticleIdxs=subParticleIdxs;
        this.particleKinds=particleKinds;
    }
    SysStruct dup(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.SysStruct){
            return new SysStruct(fullSystem,particles,superParticle,subParticleIdxs,particleKinds);
        } else if (l > PSCopyDepthLevel.SysStruct) {
            return new SysStruct(fullSystem.dup,particles.dup,superParticle.dup,
                subParticleIdxs.dup,particleKinds.dup);
        }
        return this;
    }
    SysStruct dup(){
        return dup(PSCopyDepthLevel.SysStruct);
    }
    SysStruct deepDup(){
        return dup(PSCopyDepthLevel.All);
    }
    mixin(serializeSome("dchem.sys.SysStruct",
        `fullSystem: sub mapping to the whole system
        particles: particle indexes
        superParticle: super particle, i.e. molecule for example
        subParticleIdxs: index within the super particle
        particleKinds: particle kinds`))
}

/// represent a system of particles
class ParticleSys, CopiableOjectI,Serializable
{
    SysStruct sysStruct;
    
    DynamicsVars dynVars;
    
    this(SysStruct sysStruct,DynamicsVars dynVars){
        this.sysStruct=sysStruct;
        this.dynVars=dynVars;
    }
    
    ParticleSys dup(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.None){
            return this;
        } else if (l==PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys(sysStruct.dup(l),dynVars.dup(l));
        } else if (cast(int)l>cast(int)PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys(sysStruct.dup(l),dynVars.dup(l));
        }
    }
    
    ParticleSys dup(){
        return dup(PSCopyDepthLevel.DynProperties);
    }
    
    ParticleSys deepDup(){
        return dup(PSCopyDepthLevel.All);
    }
    
    mixin(serializeSome("dchem.sys.ParticleSys",
        `sysStruct: structure of the system (particle, particle kinds,...)
        dynVars: dynamic variables (position,cell,velocities)`));
}
                                                                    
