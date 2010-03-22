module dchem.sys.ParticleSys;
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
import blip.container.BitArray;
import blip.container.Deque;
import blip.container.BulkArray;

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

char[] withParticleSysMixin(char[]tmpl,char[] variant){
    return `
    if (`~variant~`.isA!(ParticleSys!(Real))){
        `~tmpl~`!(Real)(`~variant~`.get!(ParticleSys!(Real)));
    } else if (`~variant~`.isA!(ParticleSys!(LowP))){
        `~tmpl~`!(LowP)(`~variant~`.get!(ParticleSys!(LowP)));
    } else if (`~variant~`.isA!(ParticleSys!(HighP))){
        `~tmpl~`!(HighP)(`~variant~`.get!(ParticleSys!(HighP)));
    } else {
        assert(0,"expected a ParticleSys");
    }`;
}
/// represent a group of particles with the same kind
class ParticleKind: Serializable,CopiableObjectI{
    char[] _name;
    char[] name(){ return _name; }
    void name(char[] nName){ _name=nName; }
    char[] _potential;
    char[] potential(){ return _potential; }
    void potential(char[] nPotential){ _potential=nPotential; }
    char[] symbol;
    LevelIdx _level;
    LevelIdx level(){ return _level; }
    void level(LevelIdx l){ _level=l; }
    KindIdx pKind;
    size_t position;
    size_t degreesOfFreedom;
    size_t orientation;
    size_t subParticles;
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
        metaI.addFieldOfType!(size_t)("subParticles","number of subParticles",
            SerializationLevel.debugLevel);
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    /// just for internal use
    this(){}
    
    this(char[] pName,LevelIdx pLevel,KindIdx kindIdx,char[] potential=null,char[] symbol=null,
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
    void sysStructChangedT(T)(ParticleSys!(T) p){
        p.dynVars.posStruct.addToKindDim(pKind,position);
        p.dynVars.orientStruct.addToKindDim(pKind,orientation);
        p.dynVars.dofStruct.addToKindDim(pKind,degreesOfFreedom);
    }
    
    void sysStructChanged(Variant pV){
        mixin(withParticleSysMixin("sysStructChangedT","pV"));
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(Variant p){}
    /// cell changed
    void cellChanged(Variant p){}
    /// properties to be calculated did change
    void propertiesRequestedChanged(Variant p){}
    /// properties were allocated
    void propertiesAllocated(Variant p){}
    /// will calculate the properties
    void willCalculate(Variant p){}
    /// did calculate the properties
    void didCalculate(Variant p){}
    /// did read the properties (no need to keep them in memory)
    void didReadProperties(){}
    /// request to try to reduce memory usage
    void minimizeMemory(Variant p){}
}

/// kind of the particle that repesents the whole system
class SysKind: ParticleKind{
    override char[] name(){ return "_SYSTEM_"; }
    override char[] potential() { return ""; }
    this(){}
    
    this(LevelIdx pLevel,KindIdx kindIdx,size_t position=0,size_t orientation=0, size_t degreesOfFreedom=0){
        super("_SYSTEM_",pLevel,kindIdx,"","",position,orientation,degreesOfFreedom);
    }
    
}

/// perform an operation on the segmented arrays and cell of a dynPVector
/// assumes the existence of a boolean variable named "weak" that if true suppress the exception
/// when some of the arrays are null and other aren't
/// if cell is false does the operation only on the SegmentedArray.
/// if nonEq is true checks tat the first variable (namesLocal[0].*) is different from all the others
char[] dynPVectorOp(char[][]namesLocal,char[] op,bool cell=true,bool nonEq=false)
{
    char[] res;
    char[][] els=[cast(char[])".pos",".orient",".dof"];
    if (cell) els~=".cell";
    foreach (pp;els){
        res=`
        if (`;
        foreach (i,n; namesLocal){
            if (i!=0) res~="||";
            res~=n;
            res~=" is null ";
        }
        res~=`) {`;
        if (namesLocal.length>0){
            res~=`
            if ((!weak) && (`;
            foreach (i,n; namesLocal){
                if (i!=0) res~="||";
                res~=n;
                res~=" !is null ";
            }
            res~=`)) {
                throw new Exception("non equivalent DynPVector in `~pp~`",__FILE__,__LINE__);
            }
        }`;
        } else {
            res~=`
        }`;
        }
        if (nonEq && namesLocal.length>1){
            res~=` else if (`;
            foreach (i,n;namesLocal[1..$]){
                if (i!=0) res~=" || ";
                res~=namesLocal[0]~` is `~n;
            }
            res~=`) {
            throw new Exception("equal equivalent DynPVector in `~pp~`",__FILE__,__LINE__);
        }`;
        }
        res~=` else {
            void doOp(`;
        foreach(i,n;namesLocal){
            if (i!=0) res~=",";
            res~="typeof("~n~pp~") "~n;
        }
        res~="){\n";
        res~=op;
        res~=`
            }
            doOp(`;
        foreach(i,n;namesLocal){
            if (i!=0) res~=",";
            res~=n~pp;
        }
        res~=`);
        }`;
    }
    return res;
}

/// a vector representing a state,dstate or mddstate, implements vector ops
/// the basic type to represent a state is DynamicsVars, so that forces dependent on the velocity,...
/// can be represented cleanly, but for normal conservative systems DynPVector is a useful abstraction.
struct DynPVector(T){
    Cell!(T) cell;
    /// position in 3D space
    SegmentedArray!(Vector!(T, 3)) pos;
    /// orientation (quaternions)
    SegmentedArray!(Quaternion!(T)) orient;
    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(T) dof;
    
    void axpby(V)(DynPVector!(V)x,V a,T y){
        if ((cell is x.cell) && cell !is null){
            throw new Exception("identical cells in axpby",__FILE__,__LINE__);
        }
        auto y=this;
        mixin(dynPVectorOp(["x","y"],"y.axpby(x,a,b);",true,false));
    }
    
    void opMulAssign()(T scale){
        auto x=this;
        mixin(dynPVectorOp(["x"],"x*=scale;",true,false));
    }
    void opMulAssign(V)(DynPVector!(V) y){
        auto x=this;
        mixin(dynPVectorOp(["x","y"],"x*=y;",true,false));
    }
    
    void opSliceAssign(V)(DynPVector!(V) b){
        auto a=this;
        mixin(dynPVectorOp(["a","b"],"a[]=b;",true,true));
    }
    void opSliceAssign()(T val){
        if (cell!is null)   cell[]=val;
        if (pos!is null)    pos[] =Vector!(T,3)(val);
        if (orient!is null) orient[]=Quaternion!(T).identity;
        if (dof!is null)    dof[]=val;
    }
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars(T){
    alias T dtype;
    Real potentialEnergy; /// potential energy of the system (NAN if unknown)

    // index
    // cell
    Cell!(T) cell;    /// the current cell
    Cell!(T) dcell;   /// normally not used
    Cell!(T) mddcell; /// stress tensor (often not propagated)

    /// structure of all position based arrays
    SegmentedArrayStruct posStruct;
    
    // position in 3D space
    SegmentedArray!(Vector!(T, 3)) pos;
    SegmentedArray!(Vector!(T, 3)) dpos;
    SegmentedArray!(Vector!(T, 3)) mddpos;

    /// structure of all orientation based arrays
    SegmentedArrayStruct orientStruct;

    // orientation (quaternions)
    SegmentedArray!(Quaternion!(T)) orient;
    SegmentedArray!(Quaternion!(T)) dorient;
    SegmentedArray!(Quaternion!(T)) mddorient;

    /// structure of all dof (degrees of freedom) arrays
    SegmentedArrayStruct dofStruct;

    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(T) dof;
    SegmentedArray!(T) ddof;
    SegmentedArray!(T) mdddof;
    
    /// reallocates the segmented array structures, invalidates everything
    void reallocStructs(SysStruct sys){
        posStruct=new SegmentedArrayStruct("posStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd));
        orientStruct=new SegmentedArrayStruct("orientStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd));
        dofStruct=new SegmentedArrayStruct("dofStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd));
        
           pos=null;
          dpos=null;
        mddpos=null;
           orient=null;
          dorient=null;
        mddorient=null;
           dof=null;
          ddof=null;
        mdddof=null;
    }
    
    enum DynPVCreation{
        NormalPtrCopy,
        DuplicateCell,
        NullCell,
    }
    /// position vector
    DynPVector!(T) x(DynPVCreation how=DynPVCreation.NormalPtrCopy){
        DynPVector!(T) res;
        switch(how){
        case DynPVCreation.NormalPtrCopy:
            res.cell=cell;
            break;
        case DynPVCreation.DuplicateCell:
            res.cell=cell.dup();
            break;
        case DynPVCreation.NullCell:
            res.cell=null;
        }
        res.pos=pos;
        res.orient=orient;
        res.dof=dof;
        return res;
    }
    /// velocities vector
    DynPVector!(T) dx(DynPVCreation how=DynPVCreation.NormalPtrCopy){
        DynPVector!(T) res;
        switch(how){
        case DynPVCreation.NormalPtrCopy:
            res.cell=dcell;
            break;
        case DynPVCreation.DuplicateCell:
            res.cell=dcell.dup();
            break;
        case DynPVCreation.NullCell:
            res.cell=null;
        }
        res.pos=dpos;
        res.orient=dorient;
        res.dof=ddof;
        return res;
    }
    /// forces vector
    DynPVector!(T) mddx(DynPVCreation how=DynPVCreation.NormalPtrCopy){
        DynPVector!(T) res;
        switch(how){
        case DynPVCreation.NormalPtrCopy:
            res.cell=mddcell;
            break;
        case DynPVCreation.DuplicateCell:
            res.cell=mddcell.dup();
            break;
        case DynPVCreation.NullCell:
            res.cell=null;
        }
        res.pos=mddpos;
        res.orient=mddorient;
        res.dof=mdddof;
        return res;
    }
    /// ensures that the positions are allocated
    void checkX(){
        assert(posStruct!is null);
        assert(orientStruct!is null);
        assert(dofStruct!is null);
        if (pos is null)    pos=new SegmentedArray!(Vector!(T,3))(posStruct);
        if (orient is null) orient=new SegmentedArray!(Quaternion!(T))(orientStruct);
        if (dof is null)    dof=new SegmentedArray!(T)(dofStruct);
    }
    /// ensures that the velocities are allocated
    void checkDx(){
        assert(posStruct!is null);
        assert(orientStruct!is null);
        assert(dofStruct!is null);
        if (dpos is null)    dpos=new SegmentedArray!(Vector!(T,3))(posStruct);
        if (dorient is null) dorient=new SegmentedArray!(Quaternion!(T))(orientStruct);
        if (ddof is null)    ddof=new SegmentedArray!(T)(dofStruct);
    }
    /// ensures that the forces are allocated
    void checkMddx(){
        assert(posStruct!is null);
        assert(orientStruct!is null);
        assert(dofStruct!is null);
        if (mddpos is null)    mddpos=new SegmentedArray!(Vector!(T,3))(posStruct);
        if (mddorient is null) mddorient=new SegmentedArray!(Quaternion!(T))(orientStruct);
        if (mdddof is null)    mdddof=new SegmentedArray!(T)(dofStruct);
    }
    
    void opSliceAssign(V)(ref DynamicsVars!(V) d2){
        potentialEnergy=d2.potentialEnergy;
        if (this.cell!is null && d2.cell !is null)
            this.cell[]=d2.cell;
        if (this.dcell!is null && d2.dcell !is null)
            this.dcell[]=d2.dcell;
        if (this.mddcell!is null && d2.mddcell !is null)
            this.mddcell[]=d2.mddcell;
        if (this.pos!is null && d2.pos !is null)
            d2.pos.dupTo(pos);
        if (this.dpos!is null && d2.dpos !is null)
            d2.dpos.dupTo(dpos);
        if (this.mddpos!is null && d2.mddpos !is null)
            d2.mddpos.dupTo(mddpos);
        if (this.orient!is null && d2.orient !is null)
            d2.orient.dupTo(orient);
        if (this.dorient!is null && d2.dorient !is null)
            d2.dorient.dupTo(dorient);
        if (this.mddorient!is null && d2.mddorient !is null)
            d2.mddorient.dupTo(mddorient);
        if (this.dof!is null && d2.dof !is null)
            d2.dof.dupTo(dof);
        if (this.ddof!is null && d2.ddof !is null)
            d2.ddof.dupTo(ddof);
        if (this.mdddof!is null && d2.mdddof !is null)
            d2.mdddof.dupTo(mdddof);
    }
    void opSliceAssign()(T val){
        if (cell!is null)       cell[]=val;
        if (dcell!is null)      dcell[]=val;
        if (mddcell!is null)    mddcell[]=val;
        if (pos!is null)        pos[]=val;
        if (dpos!is null)       dpos[]=val;
        if (mddpos!is null)     mddpos[]=val;
        if (orient!is null)     orient[]=val;
        if (dorient!is null)    dorient[]=val;
        if (mddorient!is null)  mddorient[]=val;
        if (dof!is null)        dof[]=val;
        if (ddof!is null)       ddof[]=val;
        if (mdddof!is null)     mdddof[]=val;
    }
    void axpby(V)(DynamicsVars!(V) v,V a,T b){
        x.axpby!(V)(v.x,a,b);
        dx.axpby!(V)(v.dx,a,b);
        mddx.axpby!(V)(v.mddx,a,b);
    }
    void opMulAssign(V)(DynamicsVars!(V) v){
        x*=v.x;
        dx*=v.dx;
        mddx*=v.mddx;
    }
    void opMulAssign()(T v){
        x*=v;
        dx*=dv;
        mddx*=mddv;
    }
    DynamicsVars!(V) dupT(V=T)(PSCopyDepthLevel level){
        if (level>=PSCopyDepthLevel.DynProperties){
            DynamicsVars!(V) res;
            res.potentialEnergy=potentialEnergy;
            res.cell=cell.dup!(V);
            res.dcell=dcell.dup!(V);
            res.mddcell=mddcell.dup!(V);
            res.pos=pos.dupT!(Vector!(V, 3))();
            res.dpos=dpos.dupT!(Vector!(V, 3))();
            res.mddpos=mddpos.dupT!(Vector!(V, 3))();
            res.orient=orient.dupT!(Quaternion!(V))();
            res.dorient=orient.dupT!(Quaternion!(V))();
            res.mddorient=orient.dupT!(Quaternion!(V))();
            res.dof=dof.dupT!(V)();
            res.ddof=ddof.dupT!(V)();
            res.mdddof=mdddof.dupT!(V)();
            return res;
        }
        return *this;
    }
    DynamicsVars dup(PSCopyDepthLevel level){
        return dupT!(T)(level);
    }
    
    mixin(serializeSome("dchem.sys.DynamicsVars("~T.stringof~")","potentialEnergy|cell|pos|dpos|mddpos|orient|dorient|mddorient|dof|ddof|mdddof"));
    mixin printOut!();
}

/// represent the structure of a system of particles
class SysStruct: CopiableObjectI,Serializable
{
    char[] name; /// name of the structure
    SubMapping fullSystem; /// LocalIndex is this system struct
    SubMapping externalOrder; /// LocalIndex is the external one
    KindRange[] levels; /// disjoint KindRange at each level
    SegmentedArrayStruct particlesStruct;
    SegmentedArray!(PIndex) particles; // make it non explicit? it would spare quite some memory...
    SegmentedArray!(PIndex) superParticle;
    SegmentedArrayStruct subParticlesStruct;
    SegmentedArray!(PIndex) subParticles;
    SegmentedArrayStruct particleKindsStruct;
    SegmentedArray!(ParticleKind) particleKinds;
    
    this(char[] name,SubMapping fullSystem,SubMapping externalOrder,KindRange[] levels,
        SegmentedArrayStruct particlesStruct,SegmentedArray!(PIndex) particles,
        SegmentedArray!(PIndex) superParticle,SegmentedArrayStruct subParticlesStruct,
        SegmentedArray!(PIndex) subParticles,
        SegmentedArrayStruct particleKindsStruct,SegmentedArray!(ParticleKind) particleKinds)
    {
        this.name=name;
        this.fullSystem=fullSystem;
        this.externalOrder=externalOrder;
        this.levels=levels;
        this.particlesStruct=particlesStruct;
        this.particles=particles;
        this.superParticle=superParticle;
        this.subParticlesStruct=subParticlesStruct;
        this.subParticles=subParticles;
        this.particleKindsStruct=particleKindsStruct;
        this.particleKinds=particleKinds;
    }
    this(){ }
    /// kinds as a simple BulkArray (easier indexing)
    BulkArray!(ParticleKind) kinds(){
        return particleKinds.data();
    }
    SysStruct dup(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.SysStruct){
            return new SysStruct(name, fullSystem, externalOrder, levels, particlesStruct, particles,
                superParticle, subParticlesStruct, subParticles, particleKindsStruct, particleKinds);
        } else if (l > PSCopyDepthLevel.SysStruct) {
            return new SysStruct(name, fullSystem, externalOrder, levels.dup, particlesStruct, particles.dup,
                superParticle.dup, subParticlesStruct, subParticles.dup, particleKindsStruct, particleKinds.dup);
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
        `name: the name of the system
        fullSystem: sub mapping to the whole system
        externalOrder: order of the external files
        levels: kind ranges of the various levels
        particles: particle indexes
        superParticle: super particle, i.e. molecule for example
        subParticles: the subparticles of each particle
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
class ParticleSys(T): CopiableObjectI,Serializable
{
    alias T dtype;
    char[] name; /// name of the particle system
    ulong iteration;
    NotificationCenter nCenter;
    SysStruct sysStruct;
    
    DynamicsVars!(T) dynVars;
    
    /// internal use
    this(){}
    /// constructor
    this(ulong iter,char[] name,SysStruct sysStruct,NotificationCenter nCenter,DynamicsVars!(T) dynVars){
        this.iteration=iter;
        this.name=name;
        this.sysStruct=sysStruct;
        this.dynVars=dynVars;
        this.nCenter=nCenter;
    }
    this(ulong iter,char[] name,SysStruct sysStruct,NotificationCenter nCenter){
        this.iteration=iter;
        this.name=name;
        this.sysStruct=sysStruct;
        this.nCenter=nCenter;
    }
    
    void zero(bool x=true,bool dx=true,bool mddx=true){
        if (x) dynVars.x[]=0;
        if (dx) dynVars.dx[]=0;
        if (mddx) dynVars.mddx[]=0;
    }
    
    ParticleSys!(V) dupT(V=T)(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.None){
            return this;
        } else if (l==PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys!(V)(iteration,name,sysStruct.dup(l),null,dynVars.dupT!(V)(l));
        } else if (cast(int)l>cast(int)PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys!(V)(iteration,name,sysStruct.dup(l),null,dynVars.dupT!(V)(l));
        }
    }
    
    ParticleSys dup(PSCopyDepthLevel l){
        return dupT!(T)(l);
    }
    
    ParticleSys dup(){
        return dupT!(T)(PSCopyDepthLevel.DynProperties);
    }
    
    ParticleSys deepdup(){
        return dup(PSCopyDepthLevel.All);
    }
    
    void reallocStructs(){
        dynVars.reallocStructs(sysStruct);
    }
    /// ensures that the positions are allocated
    void checkX(){
        dynVars.checkX();
    }
    /// ensures that the velocities are allocated
    void checkDx(){
        dynVars.checkDx();
    }
    /// ensures that the forces are allocated
    void checkMddx(){
        dynVars.checkMddx();
    }
    
    /// system structure changed (particle added/removed, kinds added/removed)
    /// the segmented array structs should be initialized, and modifiable.
    /// positions,... are not yet valid
    void sysStructChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.sysStructChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("sysStructChanged",Variant(this));
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.positionsChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("positionsChanged",Variant(this));
    }
    /// cell changed
    void cellChanged(){
        foreach(pKind;sysStruct.particleKinds.pLoop){
            pKind.cellChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("cellChanged",Variant(this));
    }
    /// copy op
    void opSliceAssign(V)(ParticleSys!(V) p){
        iteration=p.iteration;
        sysStruct=p.sysStruct;
        dynVars[]=p.dynVars;
    }
    
    mixin(serializeSome("dchem.sys.ParticleSys",
        `sysStruct: structure of the system (particle, particle kinds,...)
        dynVars: dynamic variables (position,cell,velocities)`));
    mixin printOut!();
}

/// keeps a history of the previous steps (particleSys)
class HistoryManager(T){
    size_t nHistory; /// number of steps to keep
    BitArray keepH; /// things to keep in history
    enum bPos{ /// bit positions in keepH
        cell=0,
        pos=1,
        dpos=2,
        mddpos=3,
        orient=4,
        dorient=5,
        mddorient=6,
        dof=7,
        ddof=8,
        mdddof=9
    }
    Deque!(Real) historyE;
    Deque!(ParticleSys!(T)) history; /// place to keep the history
    this(size_t nHistory=1,BitArray keepH=BitArray([true,true,false,false,true,false,false,true,false,false])){
        this.nHistory=nHistory;
        this.keepH=keepH;
        history=new Deque!(ParticleSys!(T))(nHistory);
    }
    void addToHistory(V)(ParticleSys!(V) p,Real e){
        if (nHistory==0) return;
        if (history.length<nHistory){
            DynamicsVars!(T) dVars;
            if (keepH[bPos.cell])
                dVars.cell=p.cell.dup;
            if (keepH[bPos.pos])
                dVars.pos=p.pos.dupT!(Vector!(T,3));
            if (keepH[bPos.dpos])
                dVars.dpos=p.dpos.dupT!(Vector!(T,3));
            if (keepH[bPos.mddpos])
                dVars.mddpos=p.mddpos.dupT!(Vector!(T,3));
            if (keepH[bPos.orient])
                dVars.orient=p.orient.dupT!(Quaternion!(T));
            if (keepH[bPos.dorient])
                dVars.dorient=p.dorient.dupT!(Quaternion!(T));
            if (keepH[bPos.mddorient])
                dVars.mddorient=p.mddorient.dupT!(Quaternion!(T));
            if (keepH[bPos.dof])
                dVars.dof=p.dof.dupT!(T);
            if (keepH[bPos.ddof])
                dVars.ddof=p.ddof.dupT!(T);
            if (keepH[bPos.mdddof])
                dVars.mdddof=p.mdddof.dupT!(T);
            ParticleSys!(T) n=new ParticleSys!(T)(collectAppender(delegate(CharSink s){
                s("history-"); writeOut(history.length); }),
                p.sysStruct,dVars,cast(NotificationCenter)null);
            history.pushFront(n);
            historyE.pushFront(e);
        } else {
            auto nO=history.popEnd();
            nO[]=p;
            history.pushFront(nO);
            historyE.popFront(e);
            historyE.pushFront(e);
        }
    }
}
