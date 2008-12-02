/++
 +  the simulation cell
 +/
module sysStruct.cell;
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
    Real[3][3] h,hInv;
    Real[3] x0;
}

/++
 + Particle
 + minimal information on a particle, this is a structure, and should be kept as small 
 + as possible
 + - ikind: index of the kind of the current particle
 + - imol_struct: index of the current molecule (-1 if none)
 + - subPIdx: index of the current particle in its superparticle
 +/
struct Particle
{
    int ikind;
    int imol;
    int subPIdx;
}

/++
 + ParticleKind
 + class that describes the kind of a particle. As only a few kinds per
 + simulation are expected this class can be derived, and can contain more information.
 + The information here should still be as much as possible independent of the
 + configuration/system, so that it can be shared by two systems with different number
 + of particles.
 + A particle can also be composite (molecule or fragment).
 + - name: name of the current kind, should be unique
 + - kind_type: can be used to quickly find out of which subclass of kind the actual
 +   object is instance
 + - level: level of the current particle (0-> atom 1-> molcule of atoms... (whole system))
 + - mass: mass of the current particle in atomic units
 + - charge: charge of the actual particle (e-units)
 + - subParticleKinds: sub particles (if this is a composite particle)
 + - dynProp: dynamic properties of the actual particle (position, speed, orientation,...)
 + - constraints: the constraints on the elements of this particle
 +/
class ParticleKind
{
    char[] name;
    int kind_type;
    int level;
    ParticleKind[] subParticleKinds;
    DynPPropKind[] dynPPropKinds;
    Constraints[] constaints;
    bool isAtom(){
        return subParticleKinds is null || subParticleKinds.length==0;
    }
    this(char[] name, int kind_type, int level, ParticleKind[] subParticleKinds,
        DynPPropKind[] dynPPropKinds, Constraints[] constaints)
    {
        this.name=name.dup;
        this.kind_type=kind_type;
        this.level=level;
        this.subParticleKinds=subParticleKinds;
        this.dynPPropKinds=dynPPropKinds;
        foreach(prop;dynPPropKinds){
            if (prop.name=="v") v_=prop;
            if (prop.name=="pos") pos_=prop;
            if (prop.name=="mass") mass_=prop;
            if (prop.name=="charge") charge_=prop;
        }
        this.constaints=constaints;
    }
    DynPPropKind v(){ return v_; }
    DynPPropKind pos() { return pos_; }
    DynPPropKind mass() { return mass_; }
    DynPPropKind charge() { return charge_; }
private:
    DynPPropKind v_;
    DynPPropKind pos_;
    DynPPropKind mass_;
    DynPPropKind charge_;
}

interface LocalProperty(T){
    T opIndex(idxType i);
    void opIndexAssign(idxType i,T val);
}

interface GlobalProperty(T){
    T value();
    void value(T val);
}

/++
+ Constraint is the base class that describes the constraints of a system
+/
class Constraint
{
    int[] involvedBottomAtoms();
}

/++
+ storage for the property
+/
enum StorageType{
    DENSE_REAL;
    DENSE_REAL1D;
    DENSE_REAL2D;
    DENSE_REAL3D;
}

/+
+ DynPPropKind is the base class to describe the dynamic properties of a particle
+/
class DynPPropKind
{
    char[] name;
    int prop_id;
    bool global;
    StorageType storage; /// relative to the single particle, even if global==false
    Variant defaultValue;
    DynamicPProperty allocStorage(int nparticles);
}

/+
+ DynSysPropKind is the base class to describe the dynamic properties of the system
+/
class DynSysPropKind
{
    char[] name;
    int prop_id;
    StorageType storage;
    Variant defaultValue;
    DynamicSysProperty allocStorage(ParticleSys pSys);
}

/++
+ DynPProperty stores the value of a dynamic property
+/
class DynPProperty
{
    DynPPropKind kind;
    Variant storage();
    void copyFrom(DynPProperty p);
    void resetToDefault();
    this(DynPPropKind kind){
        this.kind=kind;
    }
}

/++
+ DynRProperty stores a particle property of type T (global or not)
+/
class DynGenProperty(T): DynPProperty
{
    T valueP(uint iParticle);
    void value(T newVal);
}

/++
+ DynGlobalProperty stores a global (of ParticleKind) property of type T
+/
class DynGlobalProperty(T): DynGenProperty(T)
{
    this(DynPPropKind kind){ this(kind,kind.defaultValue.get!(T)());}
    this(DynPPropKind kind,T initialVal){
        super(kind);
        value_=initialVal;
    }
    Variant storage(){ return Variant(value_); }
    T value(){
        return val_;
    }
    override T valueP(uint iParticle){
        return value_;
    }
    override void value(T newVal){
        value_=newVal;
    }
private:
    T value_;
}

/++
+ DynLocalProperty stores a global (same for all particles) property of type T
+/
class DynGlobalProperty(T): DynGenProperty(T)
{
    this(DynPPropKind kind){ this(kind,kind.defaultValue.get!(T)());}
    this(DynPPropKind kind,T initialVal){
        super(kind);
        value_=initialVal;
    }
    Variant storage(){ return Variant(value_); }
    override T valueP(uint iParticle){
        return value_;
    }
    override void value(T newVal){
        value_=newVal;
    }
private:
    T value_;
}

/++
+ DynLocalProperty stores a local (of Particle) property of type T
+/
class DynLocalProperty(T): DynGenProperty(T)
{
    this(DynPPropKind kind){ this(kind,kind.defaultValue.get!(T)());}
    this(DynPPropKind kind,T initialVal){
        super(kind);
        value_=initialVal;
    }
    Variant storage(){ return Variant(value_); }
    override T valueP(uint iParticle){
        return value_;
    }
    override void value(T newVal){
        value_=newVal;
    }
private:
    T value_;
}

/++
+ DynLocalNAProperty stores a local (of Particle) property of type T
+/
class DynLocalNAProperty(T,int rank): DynGenProperty(NArray!(T,rank+1))
{
    this(DynPPropKind kind){ this(kind,kind.defaultValue.get!(T)());}
    this(DynPPropKind kind,T initialVal){
        super(kind);
        value_=initialVal;
    }
    Variant storage(){ return Variant(value_); }
    override T valueP(uint iParticle){
        return value_;
    }
    override void value(T newVal){
        value_=newVal;
    }
private:
    T value_;
}

/++
+ DynGRProperty stores a global Real
+/
class DynGRProperty: DynPProperty
{
    Real value_;
    this(Real defaultVal){
        super(Variant(defaultVal));
    }
    Real value(uint iParticle){
        return value_;
    }
    Real value(){
        return value_;
    }
    void value(Real newVal){
        value_=newVal;
    }
}

/++
 + ParticleSys describes a system of particles
 +/
struct ParticleSys
{
    ParticleKinds[][] pKinds;
    Particle[][] particles;
    ParticleLevel[] pLevels;
    DynSysPropKind[] sysProperties;
    DynamicPProperty[][][] particleConf;
    DynamicSysProperty[] sysConf;
}
