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
import blip.parallel.smp.WorkManager;
import blip.t.math.Math:sqrt;

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
    
    this(char[] pName,LevelIdx pLevel,KindIdx kindIdx,char[] symbol=null,char[] potential=null,
        size_t position=1,size_t orientation=0, size_t degreesOfFreedom=0){
        this._name=pName;
        this.symbol=symbol;
        this._potential=potential;
        this._level=pLevel;
        this.pKind=kindIdx;
        this.position=position;
        this.degreesOfFreedom=degreesOfFreedom;
        this.orientation=orientation;
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
/// if nonEq is true checks that the first variable (namesLocal[0].*) is different from all the others
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

/// a pool for dynamical properties
class DynPPool(T){
    this(SegmentedArrayStruct posStruct,SegmentedArrayStruct orientStruct,SegmentedArrayStruct dofStruct){
        poolPos=new SegArrPool!(Vector!(T,3))(posStruct);
        poolOrient=new SegArrPool!(Quaternion!(T))(orientStruct);
        poolDof=new SegArrPool!(T)(dofStruct);
    }
    
    SegArrPool!(Vector!(T,3)) poolPos;
    SegArrPool!(Quaternion!(T)) poolOrient;
    SegArrPool!(T) poolDof;
    /// returns a postition like array
    SegmentedArray!(Vector!(T,3)) newPos(){
        return poolPos.getObj();
    }
    /// returns an orientation like array
    SegmentedArray!(Quaternion!(T)) newOrient(){
        return poolOrient.getObj();
    }
    /// returns a dof like array
    SegmentedArray!(T) newDof(){
        return poolDof.getObj();
    }
    void flush(){
        poolPos.flush();
        poolOrient.flush();
        poolDof.flush();
    }
    void stopCaching(){
        poolPos.stopCaching();
        poolOrient.stopCaching();
        poolDof.stopCaching();
    }
}

/// a vector representing a state,dstate or mddstate, implements vector ops
/// the basic type to represent a state is DynamicsVars, so that forces dependent on the velocity, forces,...
/// can be represented cleanly, but for normal conservative systems DynPVector is a useful abstraction.
struct DynPVector(T){
    Cell!(T) cell;
    /// position in 3D space
    SegmentedArray!(Vector!(T, 3)) pos;
    /// orientation (quaternions)
    SegmentedArray!(Quaternion!(T)) orient;
    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(T) dof;
    
    /// returns a copy with a nullified cell
    DynPVector nullCell(){
        auto res=*this;
        res.cell=null;
        return res;
    }
    /// returns a copy with a duplicated cell
    DynPVector dupCell(bool shouldThrow=true){
        auto res=*this;
        if (cell!is null){
            res.cell=cell.dup;
        } else if (shouldThrow){
            throw new Exception("cannot duplicate null cell",__FILE__,__LINE__);
        }
        return res;
    }
    DynPVector *mkNullCell(){
        cell=null;
        return this;
    }
    void clear(){
        cell=null;
        pos=null;
        orient=null;
        dof=null;
    }
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
        if (pos!is null)    pos[] =Vector!(T,3)(val,val,val);
        if (orient!is null) orient[]=Quaternion!(T).identity; // not really possible, not a vector...
        if (dof!is null)    dof[]=val;
    }
    
    T opIndex(size_t idx){
        if (pos!is null){
            auto len=3*pos.data.length;
            if (idx<len){
                return pos.data[idx/3][idx%3];
            }
            idx-=len;
        }
        if (cell!is null){
            if (idx<9){
                return cell.h.cell[idx];
            }
            idx-=9;
        }
        if (dof!is null){
            auto len=dof.data.length;
            if (idx<len){
                return dof.data[idx];
            }
            idx-=len;
        }
        if (orient!is null){
            auto len=3*orient.data.length;
            if (idx<len){
                return orient.data[idx/3].xyzw.cell[idx%3];
            }
        }
        throw new Exception("index out of bounds",__FILE__,__LINE__);
    }

    void opIndexAssign(T val,size_t idx){
        if (pos!is null){
            auto len=3*pos.data.length;
            if (idx<len){
                auto v=pos.data.ptrI(idx/3);
                v.opIndexAssign(val,idx%3);
                return;
            }
            idx-=len;
        }
        if (cell!is null){
            if (idx<9){
                cell.h.cell[idx]=val;
                return;
            }
            idx-=9;
        }
        if (dof!is null){
            auto len=dof.data.length;
            if (idx<len){
                dof.data.opIndexAssign(val,idx);
                return;
            }
            idx-=len;
        }
        if (orient!is null){
            auto len=3*orient.data.length;
            if (idx<len){
                auto p=orient.data.ptrI(idx/3);
                p.xyzw.cell[idx%3]=val;
                auto n=p.xyzw.x*p.xyzw.x+p.xyzw.y*p.xyzw.y+p.xyzw.z*p.xyzw.z;
                if (n>1){
                    p.xyzw/=sqrt(n);
                    p.xyzw.w=0;
                } else {
                    p.xyzw.w=sqrt(1-n);
                }
                return;
            }
        }
        throw new Exception("index out of bounds",__FILE__,__LINE__);
    }
    
    size_t length(){
        size_t len=0;
        if (pos!is null){
            len+=3*pos.data.length;
        }
        if (cell!is null){
            len+=9;
        }
        if (dof!is null){
            len+=dof.data.length;
        }
        if (orient!is null){
            len+=3*orient.data.length;
        }
        return len;
    }
    
    U opDot(V,U=typeof(T.init+V.init))(DynPVector!(V) v2){ // could overlap the various loops...
        U res=0;
        Exception e=null;
        void doPosLoop(){
            try{
                if (e!is null) return;
                assert(v2.pos!is null,"Different vectors in opDot");
                assert(v2.pos.arrayStruct is pos.arrayStruct,"Different array structs in opDot");
                auto d1=a2NA2((cast(T*)pos.data.ptr)[3*pos.data.length]);
                auto d2=a2NA2((cast(T*)v2.pos.data.ptr)[3*v2.pos.data.length]);
                auto rAtt=dot(d1,d2);
                atomicAdd(res,cast(U)rAtt);
            } catch (Exception eTmp){
                e=new Exception("exception in doPosLoop",__FILE__,__LINE__,eTmp);
            }
        }
        void doDofLoop(){
            try{
                if (e!is null) return;
                assert(v2.dof!is null,"Different vectors in opDot");
                assert(v2.dof.arrayStruct is dof.arrayStruct,"Different array structs in opDot");
                if (dof.length==0) return;
                auto d1=a2NA2(dof.data.data);
                auto d2=a2NA2(v2.dof.data);
                auto rAtt=dot(d1,d2);
                atomicAdd(res,cast(U)rAtt);
            } catch(Exception eTmp){
                e=new Exception("exception in doDofLoop",__FILE__,__LINE__,eTmp);
            }
        }
        void doOrientLoop(){
            try{
                if (e!is null) return;
                assert(v2.orient!is null,"Different vectors in opDot");
                assert(v2.orient.arrayStruct is orient.arrayStruct,"Different array structs in opDot");
                if (orient.length==0) return;
                auto d1=NArray!(T,2)([cast(index_type)4*T.sizeof,T.sizeof],[cast(index_type)orient.data.length,3],
                    cast(index_type)0,(cast(T*)orient.data.ptr)[4*orient.data.length],0);
                auto d2=NArray!(T,2)([cast(index_type)4*T.sizeof,T.sizeof],[cast(index_type)v2.orient.data.length,3],
                    cast(index_type)0,(cast(T*)v2.orient.data.ptr)[4*v2.orient.data.length],0);
                auto rAtt=dot(d1,d2);
                atomicAdd(res,cast(U)rAtt);
            } catch(Exception eTmp){
                e=new Exception("exception in doOrientLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        Task("DynPVectorDot",
            delegate void(){
                if (pos!is null){
                    Task("DynPVectorDotPos",&doPosLoop).autorelease.submitYield();
                } else {
                    assert(v2.pos is null,"different vectors in dot");
                }
                if (cell!is null){
                    U resTmp=0;
                    auto c1=cell.h.cell[];
                    auto c2=v2.cell.h.cell[];
                    for(int i=0;i<9;++i){
                        resTmp+=c1[i]*c2[i];
                    }
                    atomicAdd(res,resTmp);
                } else {
                    assert(v2.cell is null,"different vectors in opDot");
                }
                if (dof!is null){
                    Task("DynPVectorDotDof",&doDofLoop).autorelease.submitYield();
                } else {
                    assert(v2.dof is null,"different vectors in opDot");
                }
                if (orient!is null){
                    Task("DynPVectorDotOrient",&doOrientLoop).autorelease.submit();
                } else {
                    assert(v2.orient is null,"different vectors in opDot");
                }
            }
        ).autorelease.executeNow();
        if (e!is null) throw e;
    }
    
    int opApply(int delegate(ref T)loopBody){
        int res=0;
        Exception e=null;
        void doPosLoop(){
            if (res!=0) return;
            try{
                auto resTmp=pos.pLoop.opApply(delegate int(ref Vector!(T,3) v){ // could be flattened using a NArray on pos.data like the dot...
                    if (auto res=loopBody(v.x)) return res;
                    if (auto res=loopBody(v.y)) return res;
                    if (auto res=loopBody(v.z)) return res;
                    return 0;
                });
            } catch (Exception eTmp){
                e=new Exception("exception in doPosLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doDofLoop(){
            try{
                auto resTmp=dof.data.opApply(loopBody);
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doDofLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doOrientLoop(){
            try{
                auto resTmp=orient.data.opApply(delegate int(ref Quaternion!(T) q){ // could be flattened using a NArray on pos.data...
                    if(auto r=loopBody(q.xyzw.x)) return r;
                    if(auto r=loopBody(q.xyzw.y)) return r;
                    if(auto r=loopBody(q.xyzw.z)) return r;
                });
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doOrientLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        Task("DynPVectorOpApply",
            delegate void(){
                if (pos!is null){
                    Task("DynPVectorOpApplyPos",&doPosLoop).autorelease.submitYield();
                }
                if (cell!is null){
                    auto resTmp=cell.h.opApply(loopBody);
                    if (resTmp!=0) res=resTmp;
                }
                if (dof!is null){
                    Task("DynPVectorOpApplyDof",&doDofLoop).autorelease.submitYield();
                }
                if (orient!is null){
                    Task("DynPVectorOpApplyOrient",&doOrientLoop).autorelease.submit();
                }
            }
        ).autorelease.executeNow();
        if (e) throw e;
        return res;
    }
    
    DynPVector emptyCopy(){
        DynPVector res;
        if (cell!is null) res.cell=cell.dup();
        if (pos!is null) res.pos=pos.emptyCopy();
        if (orient!is null) res.orient=orient.emptyCopy();
        if (dof!is null) res.dof=dof.emptyCopy();
        return res;
    }
    void giveBack(){
        cell=null;
        if (pos!is null) pos.giveBack();
        pos=null;
        if (orient!is null) orient.giveBack();
        orient=null;
        if (dof!is null) dof.giveBack();
        dof=null;
    }
    DynPVector dup(){
        DynPVector res;
        if (cell!is null) res.cell=cell.dup();
        if (pos!is null) res.pos=pos.dup();
        if (orient!is null) res.orient=orient.dup();
        if (dof!is null) res.dof=dof.dup();
        return res;
    }
    DynPVector!(V) dupT(V)(){
        DynPVector!(V) res;
        if (cell!is null) res.cell=cell.dupT!(V);
        if (pos!is null) res.pos=pos.dupT!(Vector!(V,3));
        if (orient!is null) res.orient=orient.dupT!(Quaternion!(V));
        if (dof!is null) res.dof=dof.dupT!(V);
        return res;
    }
    mixin(serializeSome("dchem.sys.DynamicsVars("~T.stringof~")","cell|pos|orient|dof"));
    mixin printOut!();
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars(T){
    DynPPool!(T) pool;
    alias T dtype;
    Real potentialEnergy; /// potential energy of the system (NAN if unknown)

    /// structure of all position based arrays
    SegmentedArrayStruct posStruct;
    /// structure of all orientation based arrays
    SegmentedArrayStruct orientStruct;
    /// structure of all dof (degrees of freedom) arrays
    SegmentedArrayStruct dofStruct;
    
    /// position vector
    DynPVector!(T) x;
    /// velocities vector
    DynPVector!(T) dx;
    /// forces vector
    DynPVector!(T) mddx;
    
    /// reallocates the segmented array structures, invalidates everything
    void reallocStructs(SysStruct sys){
        posStruct=new SegmentedArrayStruct("posStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            [],SegmentedArrayStruct.Flags.None);
        orientStruct=new SegmentedArrayStruct("orientStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            [],SegmentedArrayStruct.Flags.None);
        dofStruct=new SegmentedArrayStruct("dofStruct",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            [],SegmentedArrayStruct.Flags.None);
        if (pool!is null){
            pool.flush(); // call stopCaching???
        }
        pool=new DynPPool!(T)(posStruct,orientStruct,dofStruct);
        x.clear();
        dx.clear();
        mddx.clear();
    }
    
    /// checks that the given vector is allocated
    void checkAllocDynPVect(V)(DynPVector!(V)* v){
        assert(posStruct!is null);
        assert(orientStruct!is null);
        assert(dofStruct!is null);
        static if (is(V==T)){
            if (v.pos is null)    v.pos=pool.newPos();
            if (v.orient is null) v.orient=pool.newOrient();
            if (v.dof is null)    v.dof=pool.newDof();
        } else {
            if (v.pos is null)    v.pos=new SegmentedArray!(Vector!(V,3))(posStruct);
            if (v.orient is null) v.orient=new SegmentedArray!(Quaternion!(V))(orientStruct);
            if (v.dof is null)    v.dof=new SegmentedArray!(V)(dofStruct);
        }
    }
    
    /// ensures that the positions are allocated
    void checkX(){
        checkAllocDynPVect(&x);
    }
    /// ensures that the velocities are allocated
    void checkDx(){
        checkAllocDynPVect(&dx);
    }
    /// ensures that the forces are allocated
    void checkMddx(){
        checkAllocDynPVect(&mddx);
    }
    
    void opSliceAssign(V)(ref DynamicsVars!(V) d2){
        potentialEnergy=d2.potentialEnergy;
        x[]=d2.x;
        dx[]=d2.dx;
        mddx[]=d2.mddx;
    }
    void opSliceAssign()(T val){
        potentialEnergy=0;
        x[]=val;
        dx[]=val;
        mddx[]=val;
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
            res.x=x.dupT!(V);
            res.dx=dx.dupT!(V);
            res.mddx=mddx.dupT!(V);
            return res;
        }
        return *this;
    }
    DynamicsVars dup(PSCopyDepthLevel level){
        return dupT!(T)(level);
    }
    T opDot(V)(ref DynamicsVars!(V) b){ // could be parallel...
        auto r1=x.opDot(b.x);
        auto r2=dx.opDot(b.dx);
        auto r3=mddx.opDot(b.mddx);
        return r1+r2+r3;
    }
    size_t length(){
        return x.length+dx.length+mddx.length;
    }
    /// flat indexing of the contents
    T opIndex(size_t i){
        auto len=x.length;
        if (i<len){
            return x[i];
        }
        i-=len;
        len=dx.length;
        if (i<len){
            return dx[i];
        }
        i-=len;
        assert(i<mddx.length,"index out of bounds");
        return mddx[i];
    }
    /// flat indexing of the contents
    T opIndexAssign(T val,size_t i){
        auto len=x.length;
        if (i<len){
            return x[i];
        }
        i-=len;
        len=dx.length;
        if (i<len){
            return dx[i];
        }
        i-=len;
        assert(i<mddx.length,"index out of bounds");
        return mddx[i];
    }
    mixin(serializeSome("dchem.sys.DynamicsVars("~T.stringof~")","potentialEnergy|x|dx|mddx"));
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

/// interface for hidden variables
interface HiddenVars:Serializable{
    void axpby(HiddenVars x,Real a,Real b);
    HiddenVars dup();
    HiddenVars emptyCopy();
    void opSliceAssign(Real r);
    void opSliceAssign(HiddenVars b);
    void opMulAssign(Real b);
    void opMulAssign(HiddenVars b);
    Real opDot(HiddenVars b);
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
    
    DynamicsVars!(T) dynVars; /// dynamics variables
    HiddenVars hVars; /// possibly some hidden degrees of freedom
    
    /// internal use
    this(){}
    /// constructor
    this(ulong iter,char[] name,SysStruct sysStruct,NotificationCenter nCenter,DynamicsVars!(T) dynVars,
        HiddenVars hVars=null){
        this.iteration=iter;
        this.name=name;
        this.sysStruct=sysStruct;
        this.dynVars=dynVars;
        this.nCenter=nCenter;
        this.hVars=hVars;
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
        if (hVars!is null) hVars[]=0;
    }
    
    ParticleSys!(V) dupT(V=T)(PSCopyDepthLevel l){
        if (l==PSCopyDepthLevel.None){
            return this;
        } else if (l==PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys!(V)(iteration,name,sysStruct.dup(l),null,dynVars.dupT!(V)(l),
                ((hVars!is null)?hVars.dup:null));
        } else if (cast(int)l>cast(int)PSCopyDepthLevel.PSysLevel) {
            return new ParticleSys!(V)(iteration,name,sysStruct.dup(l),null,dynVars.dupT!(V)(l),
                ((hVars!is null)?hVars.dup:null));
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
    /// positions,... are not yet available
    void sysStructChanged(){
        foreach(pKind;sysStruct.particleKinds.data){
            pKind.sysStructChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("sysStructChanged",Variant(this));
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(){
        foreach(pKind;sysStruct.particleKinds.data){
            pKind.positionsChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("positionsChanged",Variant(this));
    }
    /// cell changed
    void cellChanged(){
        foreach(pKind;sysStruct.particleKinds.data){
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
        if (hVars!is null) hVars[]=p.hVars;
    }
    
    mixin(serializeSome("dchem.sys.ParticleSys",
        `sysStruct: structure of the system (particle, particle kinds,...)
        dynVars: dynamic variables (position,cell,velocities,forces,...)
        hVars: hidden degrees of freedom`));
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
