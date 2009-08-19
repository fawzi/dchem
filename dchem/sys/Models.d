module dchem.sys.Models;

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

/+
/// Meta informations, properties of a property
interface PropertyKind: CopiableOjectI,Serializable{
    enum Distribution:int{
        Replicated,
        ParticleLocal
    }
    enum DetailLevel:int{
        KindLevel,
        ParticleLevel
    }
    enum Storage:int{ // ragged (sparse) storage to add
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
    ParticleKind particleKind();
    PropertyKind dup(PSCopyDepthLevel level);
    BulkArray!(idxType)[] kinds2particles;
}

/// generic property interface
interface GenProperty:CopiableObjectI,Serializable {
    PropertyKind kind(); /// kind of this property
    Variant storage(); /// the data
    GenProperty dup(PSCopyDepthLevel level);
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
+/

final class ParticleIterFromIdx(IType,PIType):FIteratorI!(Particle*) {
    IType it;
    BulkArray!(Particle) particles
    idxType next(){
        return particles[it.next()];
    }
    bool atEnd(){
        return it.atEnd();
    }
    int opApply(int delegate(Particle*) loopBody){
        return it.opApply(delegate int(idxType idx){ return loopBody(particles[idx]); });
    }
    ParticleIterFromIdx!(PIType,PIType) parallelLoop(size_t optimalCS){
        return new ParticleIterFromIdx!(PIType,PIType)(it.parallelLoop(optimalCS));
    }
    ParticleIterFromIdx!(PIType,PIType) parallelLoop(size_t optimalCS){
        return new ParticleIterFromIdx!(PIType,PIType)(it.parallelLoop());
    }
}

/++
 +  description of the simulation cell
 +   - periodicity: an array with 1 if that dimension is periodic, 0 if it isn't
 +   - h: matrix that maps the cell [0,1]^3 to the real cell in atomic units
 +   - h_inv: inverse of h
 +   - x0: shift of the origin (the mapping reduced points -> real points is
 +     r_real=dot(h,r_red)+x0
 +/
class Cell
{
    int[3] periodic;
    NArray!(Real,2) h,hInv;
    NArray!(Real,1) x0;
    OrthoCell orthoCell
    this(NArray!(Real,2)h,int[3] periodic,NArray!(Real,2)hInv=null){
        this.h=h;
        this.periodic=periodic;
        if (hInv is null){
            this.hInv=LinAlg.inv(h);
        } else {
            this.hInv=hInv;
        }
    }
    Cell dup(){
        return new Cell(h.dup(),periodic,hInv.dup);
    }
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars{
    // cell
    Cell cell;

    // position in 3D space (used for neigh lists/hierarchical partitioning, screening)
    SegmentedArray!(float) spos;

    // position in 3D space
    SegmentedArray!(Real) pos;
    SegmentedArray!(Real) dpos;
    SegmentedArray!(Real) ddpos;

    // orientation (quaternions)
    SegmentedArray!(Real) orient;
    SegmentedArray!(Real) dorient;
    SegmentedArray!(Real) ddorient;

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
}

/// represent the structure of a system of particles
class SysStruct, CopiableOjectI //,Serializable
{
    SubMapping fullSystem;
    SegmentedArray!(PIndex) particles;
    SegmentedArray!(PIndex) superParticle;
    SegmentedArray!(size_t) subParticleIdxs;
    SegmentedArray!(PKinds) particleKinds;
    
    SysStruct dup(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.None){
            return this;
        } else if (l==PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys(cell.dup(),pLevels.dup());
        } else if (cast(int)l>cast(int)PSCopyDepthLevel.PSysLevel) {
            auto newPS=new ParticleSys(new LevelMappings[pLevels.length]);
            foreach (i,pLevel;pLevels){
                newPS.pLevels[i]=pLevel.dup(l);
            }
        }
    }
    SysStruct dup(){
        return dup(PSCopyDepthLevel.PSysLevel);
    }
    SysStruct deepDup(){
        return dup(PSCopyDepthLevel.All);
    }
}

/// represent a system of particles
class ParticleSys, CopiableOjectI //,Serializable
{
    SysStruct sysStruct;
    
    DynamicsVars dynVars;
    
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
        return dup(PSCopyDepthLevel.PSysLevel);
    }
    ParticleSys deepDup(){
        return dup(PSCopyDepthLevel.All);
    }
    
    // mixin(expose!(NewSerializationExpose)(`cell|pLevels`));
}
                                                                    
