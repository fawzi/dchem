module dchem.sys.CoreDefs;
import dchem.sys.Requests;
import dchem.sys.PIndexes;
import blip.serialization.Serialization;
import blip.serialization.StringSerialize;
import blip.serialization.SerializationMixins;
import blip.BasicModels;
import dchem.Common;
import dchem.sys.Cell;
import dchem.sys.SubMapping;
import dchem.sys.SegmentedArray;
import blip.util.NotificationCenter;
import blip.t.core.Variant;

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
class ParticleKind: CopiableObjectI,Serializable{
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
    /// just for internal use
    this(){}
    
    this(char[] pName,LevelIdx pLevel,KindIdx kindIdx,
        size_t position=1,size_t orientation=0, size_t degreesOfFreedom=0){
        _name=pName;
        _level=pLevel;
        pKind=kindIdx;
        position=1;
        degreesOfFreedom=0;
        orientation=0;
    }
    typeof(this)dup(){
        ParticleKind res=cast(ParticleKind)this.classinfo.create();
        res.copy(this);
        return res;
    }
    typeof(this)deepdup(){
        return dup();
    }
    /// copies a particleKind
    void copy(ParticleKind p){
        _name=p._name;
        _level=p._level;
        pKind=p.pKind;
        position=p.position;
        degreesOfFreedom=p.degreesOfFreedom;
        orientation=p.orientation;
    }
    
    void serial(S)(S s){
        s.field(metaI[0],_name);
        auto ui=cast(ubyte*)&_level;
        s.field(metaI[1],*ui);
        auto us=cast(ushort*)&pKind;
        s.field(metaI[2],*us);
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
    
    mixin printOut!();
    
    // callbacks
    /// system structure changed (particle added/removed, kinds added/removed)
    /// the segmented array structs should be initialized, positions,... are not yet valid
    void sysStructChanged(ParticleSys p){
        p.dynVars.posStruct.addToKindDim(pKind,position);
        p.dynVars.orientStruct.addToKindDim(pKind,orientation);
        p.dynVars.dofStruct.addToKindDim(pKind,degreesOfFreedom);
    }
    /// position of particles changed, position,... are valid
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
    void didReadProperties(){}
    /// request to try to reduce memory usage
    void minimizeMemory(ParticleSys p){}
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars{
    // index
    // cell
    Cell cell;

    /// structure of all position based arrays
    SegmentedArrayStruct posStruct;
    
    // position in 3D space
    SegmentedArray!(vec3R) pos;
    SegmentedArray!(vec3R) dpos;
    SegmentedArray!(vec3R) ddpos;

    /// structure of all orientation based arrays
    SegmentedArrayStruct orientStruct;

    // orientation (quaternions)
    SegmentedArray!(Quat) orient;
    SegmentedArray!(Quat) dorient;
    SegmentedArray!(Quat) ddorient;

    /// structure of all dof (degrees of freedom) arrays
    SegmentedArrayStruct dofStruct;

    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(Real) dof;
    SegmentedArray!(Real) ddof;
    SegmentedArray!(Real) dddof;
    
    /// structure of the constraint array
    SegmentedArrayStruct constraintsStruct;
    // constraints
    SegmentedArray!(Constraint) constraints;

    void opSliceAssign(ref DynamicsVars d2){
        if (this.cell!is null && d2.cell !is null)
            this.cell[]=d2.cell;
        if (this.pos!is null && d2.pos !is null)
            d2.pos.dupTo(pos);
        if (this.dpos!is null && d2.dpos !is null)
            d2.dpos.dupTo(dpos);
        if (this.ddpos!is null && d2.ddpos !is null)
            d2.ddpos.dupTo(ddpos);
        if (this.orient!is null && d2.orient !is null)
            d2.orient.dupTo(orient);
        if (this.dorient!is null && d2.dorient !is null)
            d2.dorient.dupTo(dorient);
        if (this.ddorient!is null && d2.ddorient !is null)
            d2.ddorient.dupTo(ddorient);
        if (this.dof!is null && d2.dof !is null)
            d2.dof.dupTo(dof);
        if (this.ddof!is null && d2.ddof !is null)
            d2.ddof.dupTo(ddof);
        if (this.dddof!is null && d2.dddof !is null)
            d2.dddof.dupTo(dddof);
    }
    DynamicsVars dup(PSCopyDepthLevel level){
        if (level>=PSCopyDepthLevel.DynProperties){
            DynamicsVars res;
            res.cell=cell.dup;
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
    
    mixin(serializeSome("dchem.sys.DynamicsVars","cell|pos|dpos|ddpos|orient|dorient|ddorient|dof|ddof|dddof"));
    mixin printOut!();
}

/// represent the structure of a system of particles
class SysStruct: CopiableObjectI,Serializable
{
    SubMapping fullSystem;
    SegmentedArray!(PIndex) particles; // make it non explicit? it would spare quite some memory...
    SegmentedArray!(PIndex) superParticle;
    SegmentedArray!(size_t) subParticleIdxs;
    SegmentedArray!(ParticleKind) particleKinds;
    
    this(){ }
    this(SubMapping fullSystem,SegmentedArray!(PIndex) particles,
        SegmentedArray!(PIndex) superParticle,SegmentedArray!(size_t) subParticleIdxs,
        SegmentedArray!(ParticleKind) particleKinds)
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
            return new SysStruct(fullSystem,particles.dup,superParticle.dup,
                subParticleIdxs.dup,particleKinds.dup);
        }
        return this;
    }
    SysStruct dup(){
        return dup(PSCopyDepthLevel.SysStruct);
    }
    SysStruct deepdup(){
        return dup(PSCopyDepthLevel.All);
    }
    mixin(serializeSome("dchem.sys.SysStruct",
        `fullSystem: sub mapping to the whole system
        particles: particle indexes
        superParticle: super particle, i.e. molecule for example
        subParticleIdxs: index within the super particle
        particleKinds: particle kinds`));
    mixin printOut!();
}

/// represent a system of particles
///
/// startup should be as follow:
/// - create valid sysStruct with valid Kinds
/// - call sysStructChanged
/// - set cell, positions,...
/// - call cellChanged
/// - call posChanged
/// 
/// later skip the calls that are not needed (i.e. if only the positions did change,
/// call just posChanged)
class ParticleSys: CopiableObjectI,Serializable
{
    char[] name; /// name of the particle system
    ulong iteration;
    NotificationCenter nCenter;
    
    SysStruct sysStruct;
    
    DynamicsVars dynVars;
    
    /// internal use
    this(){}
    /// constructor
    this(ulong iter,char[] name,SysStruct sysStruct,DynamicsVars dynVars,NotificationCenter nCenter){
        this.iter=iter;
        this.name=name;
        this.sysStruct=sysStruct;
        this.dynVars=dynVars;
        this.nCenter=nCenter;
    }
    
    ParticleSys dup(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.None){
            return this;
        } else if (l==PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys(iter,name,sysStruct.dup(l),dynVars.dup(l),new NotificationCenter());
        } else if (cast(int)l>cast(int)PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys(iter,name,sysStruct.dup(l),dynVars.dup(l),new NotificationCenter());
        }
    }
    
    ParticleSys dup(){
        return dup(PSCopyDepthLevel.DynProperties);
    }
    
    ParticleSys deepdup(){
        return dup(PSCopyDepthLevel.All);
    }
    
    /// system structure changed (particle added/removed, kinds added/removed)
    /// the segmented array structs should be initialized, positions,... are not yet valid
    void sysStructChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.sysStructChanged(this);
        }
        if (nCenter!is null)
            nCenter.notify("sysStructChanged",Variant(this));
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.positionsChanged(this);
        }
        if (nCenter!is null)
            nCenter.notify("positionsChanged",Variant(this));
    }
    /// cell changed
    void cellChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.cellChanged(this);
        }
        if (nCenter!is null)
            nCenter.notify("cellChanged",Variant(this));
    }
    /// copy op
    void opSliceAssign(ParticleSys p){
        iteration=p.iteration;
        sysStruct=p.sysStruct;
        dynVars=p.dynVars;
    }
    
    mixin(serializeSome("dchem.sys.ParticleSys",
        `sysStruct: structure of the system (particle, particle kinds,...)
        dynVars: dynamic variables (position,cell,velocities)`));
    mixin printOut!();
}
                                                                    
