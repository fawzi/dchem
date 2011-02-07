/// dynamic variables of a system of particles, vectors, matrixes, pool, and a structure to keep all of them togheter
/// 
/// to do: put DynPVectStruct in DynPVector, and support fully contiguous layouts (pos,dof,orient,cell)
/// author: Fawzi
module dchem.sys.DynVars;
import dchem.sys.PIndexes;
import blip.serialization.Serialization;
import blip.serialization.StringSerialize;
import blip.serialization.SerializationMixins;
import blip.BasicModels;
import dchem.Common;
import dchem.sys.Cell;
import dchem.sys.SubMapping;
import dchem.sys.SegmentedArray;
import blip.parallel.smp.WorkManager;
import blip.math.Math:sqrt;
import blip.core.sync.Mutex;
import blip.narray.NArray;
import blip.core.Traits: ctfe_i2a;
import blip.sync.Atomic;
import blip.util.TemplateFu: nArgs;
import blip.io.BasicIO;
import blip.container.GrowableArray;
enum{
    XType=0,
    DxType=1,
    DualDxType=2,
}

/// this is at the same time a pool and a structure for a DynPVector
/// setup as follow: init, update *Structs, possibly consolidateStructs, allocPools, possibly consolidate (normally done by DynamicsVarsStruct.freezeStructs)
class DynPVectStruct(T){
    char[] name;
    int cellPeriod; // ugly, here just to be able to alloca a cell meaningfully
    SegmentedArrayStruct posStruct;
    SegmentedArrayStruct orientStruct;
    SegmentedArrayStruct dofStruct;
    
    SegArrMemMap!(Vector!(T,3)) poolPos;
    SegArrMemMap!(Quaternion!(T)) poolOrient;
    SegArrMemMap!(T) poolDof;
    
    /// equal if it has the same structures
    override equals_t opEquals(Object o){
        auto t=cast(DynPVectStruct)o;
        if (t is null) return 0;
        return (posStruct==t.posStruct && orientStruct==t.orientStruct && dofStruct==t.dofStruct);
    }
    /// compatible if they have the same structure
    bool compatible(V)(DynPVectStruct!(V) t){
        if (t is null) return false;
        return (posStruct==t.posStruct && orientStruct==t.orientStruct && dofStruct==t.dofStruct);
    }
    mixin(descSome("dchem.DynPVectStruct!("~T.stringof~")",`name|cellPeriod|posStruct|orientStruct|dofStruct`));
    
    /// constructor, sets up the structures
    this(char[]name,SubMapping submapping,KindRange kRange,SegmentedArrayStruct.Flags f=SegmentedArrayStruct.Flags.None){
        this.name=name;
        posStruct=new SegmentedArrayStruct(name~".posStruct",submapping,kRange,[],f);
        orientStruct=new SegmentedArrayStruct(name~".orientStruct",submapping,kRange,[],f);
        dofStruct=new SegmentedArrayStruct(name~".dofStruct",submapping,kRange,[],f);
    }
    /// constructor that uses the given structures
    this(char[]name,SegmentedArrayStruct posStruct,SegmentedArrayStruct orientStruct,SegmentedArrayStruct dofStruct){
        this.name=name;
        assert(posStruct!is null,"posStruct has to be allocated");
        assert(orientStruct!is null,"orientStruct has to be allocated");
        assert(dofStruct!is null,"dofStruct has to be allocated");
        this.posStruct=posStruct;
        this.orientStruct=orientStruct;
        this.dofStruct=dofStruct;
    }
    /// constructor for internal use only
    this(){}
    /// alloc pools, reusing baseVal pools if possible
    void allocPools(DynPVectStruct!(T) baseVal=null)
    {
        posStruct.freeze();
        orientStruct.freeze();
        dofStruct.freeze();
        if (baseVal!is null && baseVal.poolPos!is null && baseVal.poolPos.arrayStruct is posStruct){
            poolPos=baseVal.poolPos;
        } else {
            poolPos=new SegArrMemMap!(Vector!(T,3))(posStruct);
        }
        if (baseVal!is null && baseVal.poolOrient!is null && baseVal.poolOrient.arrayStruct is orientStruct){
            poolOrient=baseVal.poolOrient;
        } else {
            poolOrient=new SegArrMemMap!(Quaternion!(T))(orientStruct);
        }
        if (baseVal!is null && baseVal.poolDof!is null && baseVal.poolDof.arrayStruct is dofStruct){
            poolDof=baseVal.poolDof;
        } else {
            poolDof=new SegArrMemMap!(T)(dofStruct);
        }
    }
    /// consolidates pools
    void consolidate(DynPVectStruct!(T) baseVal)
    {
        if (baseVal!is null && baseVal.poolPos!is null && baseVal.poolPos.arrayStruct is posStruct){
            poolPos=baseVal.poolPos;
        }
        if (baseVal!is null && baseVal.poolOrient!is null && baseVal.poolOrient.arrayStruct is orientStruct){
            poolOrient=baseVal.poolOrient;
        }
        if (baseVal!is null && baseVal.poolDof!is null && baseVal.poolDof.arrayStruct is dofStruct){
            poolDof=baseVal.poolDof;
        }
    }
    /// consolidates the structures that are equivalent (struct names might get lost)
    void consolidateStructs(DynPVectStruct!(T) baseVal)
    {
        assert(posStruct!is null&&orientStruct!is null && dofStruct!is null);
        if (baseVal!is null && baseVal.posStruct == posStruct){
            posStruct=baseVal.posStruct;
        }
        if (baseVal!is null && baseVal.orientStruct == orientStruct){
            orientStruct=baseVal.orientStruct;
        }
        if (baseVal!is null && baseVal.dofStruct is dofStruct){
            dofStruct=baseVal.dofStruct;
        }
    }
    /// returns a postition like array
    SegmentedArray!(Vector!(T,3)) newPos(){
        assert(poolPos!is null,"you need to call allocPools before asking elements from the pool");
        return poolPos.getObj();
    }
    /// returns an orientation like array
    SegmentedArray!(Quaternion!(T)) newOrient(){
        assert(poolOrient!is null,"you need to call allocPools before asking elements from the pool");
        return poolOrient.getObj();
    }
    /// returns a dof like array
    SegmentedArray!(T) newDof(){
        assert(poolDof!is null,"you need to call allocPools before asking elements from the pool");
        return poolDof.getObj();
    }
    /// returns a DynPVector
    DynPVector!(T,group) allocGroup(int group)(bool allocCell=false){
        return DynPVector!(T,group).allocFromPool(this,allocCell);
    }
    /// removes all cached arrays
    void flush(){
        poolPos.flush();
        poolOrient.flush();
        poolDof.flush();
    }
    /// removes cached arrays and never cache again
    void stopCaching(){
        poolPos.stopCaching();
        poolOrient.stopCaching();
        poolDof.stopCaching();
    }
    /// convex hull of the ranges of all structures
    KindRange gRange(){
        KindRange res=posStruct.kRange;
        res.convexHull(orientStruct.kRange);
        res.convexHull(dofStruct.kRange);
        return res;
    }
    /// length of a particle based loop on elements that have this structure
    size_t length(){
        return posStruct.length+orientStruct.length+dofStruct.length;
    }
    /// length of a particle based loop on elements that have this structure
    size_t dataLength(){
        return posStruct.dataLength+orientStruct.dataLength+dofStruct.dataLength;
    }
}

/// defines the possible null handlings when doing operations on segmented arrays
/// if any operator is null no operation is performed, but extra checks are performed
/// according with the requested NullHandling
enum NullHandling{
    NullConsistent,  /// either all operators are null, or none is
    IgnoreRightNull, /// the right side might have more nulls than the left
    IgnoreAllNull,   /// if any operator is null skips the operation
}

/// perform an operation on the segmented arrays and cell of a dynPVector
/// assumes the existence of a int variable named "nullHandling" that controls the handling of null 
/// if cell is false does the operation only on the SegmentedArray.
/// if nonEq is true checks that the first variable (namesLocal[0].*) is different from all the others
/// switch to storing a pointer to a DynPVectStruct?
char[] dynPVectorOp(char[][]namesLocal,char[] op,bool cell=true,bool nonEq=false)
{
    char[] res;
    char[][] els=[".pos",".orient",".dof"];
    if (cell) els=[".pos",".orient",".dof",".cell"];
    foreach (iter,pp;els){
        res~=`
        if (`;
        foreach (i,n; namesLocal){
            if (i!=0) res~="||";
            res~=n;
            res~=pp;
            res~=" is null ";
        }
        res~=`) {`;
        if (namesLocal.length>0){
            res~=`
            if ((nullHandling==NullHandling.NullConsistent && (`;
            foreach (i,n; namesLocal){
                if (i!=0) res~="||";
                res~=n;
                res~=pp;
                res~=" !is null ";
            }
            if (namesLocal.length>1){
                res~=`))||
                (nullHandling==NullHandling.IgnoreRightNull && (`;
                foreach (i,n; namesLocal[1..$]){
                    if (i!=0) res~="||";
                    res~=n;
                    res~=pp;
                    res~=" !is null ";
                }
            }
            res~=`))) {
                auto msg=collectAppender(delegate void(CharSink sink){
                    auto s=dumper(sink);
                    s("non equivalent DynPVector in `~pp~`");`;
            foreach (i,n; namesLocal){
                res~=`
                    s("`~n~`:")(`~n~`)("\n");`;
            }
            res~=`
                });
                throw new Exception(msg,__FILE__,__LINE__);
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
                res~=namesLocal[0]~pp~` is `~n~pp;
            }
            res~=`) {
            throw new Exception("equal equivalent DynPVector in `~pp~`",__FILE__,__LINE__);
        }`;
        }
        res~=` else {
            void doOp`~ctfe_i2a(iter)~`(`;
        foreach(i,n;namesLocal){
            if (i!=0) res~=",";
            res~="typeof("~n~pp~") "~n;
        }
        res~="){\n";
        res~=op;
        res~=`
            }
            doOp`~ctfe_i2a(iter)~`(`;
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
/// the basic type to represent a state is DynamicsVars, so that forces dependent on the velocity, forces,...
/// can be represented cleanly, but for normal conservative systems DynPVector is a useful abstraction.
/// group is used just to subdivide DynPVectors in typechecked incompatible groups (position (0), and derivatives (1) for example)
///
/// an important thing about DynPVector is if to store their structure or not.
/// at the moment it is nt stored, this decision might change in the future
struct DynPVector(T,int group){
    alias T dtype;
    enum{ vGroup=group }

    /// cell
    Cell!(T) cell;
    /// position in 3D space
    SegmentedArray!(Vector!(T, 3)) pos;
    /// orientation (quaternions)
    SegmentedArray!(Quaternion!(T)) orient;
    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(T) dof;
    
    /// the vector is dummy
    bool isDummy(){
        return cell is null && pos is null && orient is null && dof is null;
    }
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
    /// cast between different groups, obviously if the groups are structurally different you might get problems
    /// later on (but as different group have exactly the same content type this is still safe memorywise)
    DynPVector!(T,g2) toGroup(int g2)(){
        union GConv{
            DynPVector v1;
            DynPVector!(T,g2) v2;
        }
        GConv a;
        a.v1=*this;
        return a.v2;
    }
    void clear(){
        cell=null;
        pos=null;
        orient=null;
        dof=null;
    }
    void opBypax(V)(V x,T a=1,T b=1){
        static assert(is(V==DynPVector!(V.dtype,group)),"opBypax only between DynPVectors of the same group, not "~V.stringof);
        alias NullHandling.NullConsistent nullHandling; // relax? then lhs might need a scaling if the rhs is null
        auto y=this;
        if (cell!is null){
            if (x.cell !is null){
                bypax(y.cell,x.cell,a,b);
            } else {
                cell=cell*b;
            }
        }
        mixin(dynPVectorOp(["x","y"],"bypax(y,x,a,b);",false,false));
    }
    void opMulAssign()(T scale){
        alias NullHandling.IgnoreAllNull nullHandling;
        auto x=this;
        if (cell !is null){
            cell=cell*scale;
        }
        mixin(dynPVectorOp(["x"],"x*=scale;",false,false));
    }
    void opMulAssign(V)(DynPVector!(V,group) y){
        alias NullHandling.IgnoreRightNull nullHandling;
        auto x=this;
        if (cell !is null && y.cell !is null){
            cell=cell*y.cell;
        } else {
            assert(y.cell is null);
        }
        mixin(dynPVectorOp(["x","y"],"x*=y;",false,false));
    }
    /// outer product supported only if t or u are scalars
    static V outerOp(T,U,V,R,M)(T t,U u,V v,R scaleA,M scaleRes){
        static if(is(T.dtype)){
            v.opBypax(t,u*scaleA,scaleRes);
        } else {
            v.opBypax(u,t*scaleA,scaleRes);
        }
        return v;
    }
    void opSliceAssignT(V,int gg)(DynPVector!(V,gg) b){
        static assert(gg==group,"slice assign only within the same group");
        alias NullHandling.IgnoreRightNull nullHandling;
        auto a=this;
        mixin(dynPVectorOp(["a","b"],"a[]=b;",false,true));
        static if (is(T==V)){
            if (b.cell!is null) cell=b.cell;
        } else {
            if (b.cell!is null) cell=b.cell.dupT!(T)();
        }
    }
    void opSliceAssignEl(T val){
        if (cell!is null)   cell[]=val;
        if (pos!is null)    pos[] =Vector!(T,3)(val,val,val);
        if (orient!is null) orient[]=Quaternion!(T).identity; // not really possible, not a vector...
        if (dof!is null)    dof[]=val;
    }
    void opSliceAssign(V)(V v){
        static if (is(typeof(opSliceAssignEl(v)))){
            opSliceAssignEl(v);
        } else static if (is(typeof(this.opSliceAssignT(v)))){
            this.opSliceAssignT(v);
        } else static if (is(typeof(v.copyTo(*this)))){
            v.copyTo(*this);
        } else {
            static assert(0,"cannot assign from "~V.stringof~" to DynPVector!("~T.stringof~","~ctfe_i2a(group)~")");
        }
    }
    T opIndex(size_t idx){
        if (pos!is null){
            auto len=3*pos.dataLength;
            assert(pos.contiguous,"not implemented for non contiguous arrays");
            if (idx<len){
                return pos.support[idx/3][idx%3];
            }
            idx-=len;
        }
        if (dof!is null){
            assert(dof.contiguous,"not implemented for non contiguous arrays");
            auto len=dof.dataLength;
            if (idx<len){
                return dof.support[idx];
            }
            idx-=len;
        }
        if (orient!is null){
            assert(pos.contiguous,"not implemented for non contiguous arrays");
            auto len=3*orient.dataLength;
            if (idx<len){
                return orient.support[idx/3].xyzw.cell[idx%3];
            }
        }
        if (cell!is null){
            if (idx<9){
                return cell.h.cell[idx];
            }
            idx-=9;
        }
        throw new Exception("index out of bounds",__FILE__,__LINE__);
    }

    void opIndexAssign(T val,size_t idx){
        size_t idxOrig=idx;
        if (pos!is null){
            assert(pos.contiguous,"not implemented for non contiguous arrays");
            auto len=3*pos.dataLength;
            if (idx<len){
                auto v=pos.support.ptrI(idx/3);
                v.opIndexAssign(val,idx%3);
                return;
            }
            idx-=len;
        }
        if (dof!is null){
            assert(dof.contiguous,"not implemented for non contiguous arrays");
            auto len=dof.dataLength;
            if (idx<len){
                dof.support.opIndexAssign(val,idx);
                return;
            }
            idx-=len;
        }
        if (orient!is null){
            assert(orient.contiguous,"not implemented for non contiguous arrays");
            auto len=3*orient.dataLength;
            if (idx<len){
                auto p=orient.support.ptrI(idx/3);
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
        if (cell!is null){
            if (idx<9){
                cell=cell.dup(); // eccessive copying??? introduce ref counting and avoid copying if unique?
                cell.h.cell[idx]=val;
                cell.hChanged();
                return;
            }
            idx-=9;
        }
        throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)(idxOrig)(" index out of bounds ")(0)("..")(length);
            }),__FILE__,__LINE__);
    }
    
    size_t length(){
        size_t len=0;
        if (pos!is null){
            len+=3*pos.dataLength;
        }
        if (cell!is null){
            len+=9;
        }
        if (dof!is null){
            len+=dof.dataLength;
        }
        if (orient!is null){
            len+=3*orient.dataLength;
        }
        return len;
    }
    
    struct DotOpStruct(V,U){
        U res=0;
        Exception e=null;
        DynPVector *v1;
        V v2;
        
        void addRes(U inc){
            static if(is(typeof(atomicAdd(this.res,inc)))){
                atomicAdd(this.res,inc);
            } else {
                synchronized{
                    this.res+=inc;
                }
            }
        }
        void doPosLoop(){
            try{
                if (this.e!is null) return;
                assert(this.v2.pos!is null && this.v2.pos!is null,"Invalid vectors in opDot");
                assert(this.v2.pos.arrayStruct == this.v1.pos.arrayStruct,"Different array structs in opDot");
                assert(this.v1.pos.contiguous && this.v2.pos.contiguous,"not implemented for non contiguous arrays");
                auto d1=a2NA(this.v1.pos.support.basicData);
                auto d2=a2NA(this.v2.pos.support.basicData);
                auto rAtt=dot(d1,d2);
                addRes(cast(U)rAtt);
            } catch (Exception eTmp){
                this.e=new Exception("exception in doPosLoop",__FILE__,__LINE__,eTmp);
            }
        }
        void doDofLoop(){
            try{
                if (this.e!is null) return;
                assert(this.v2.dof!is null && this.v1.dof!is null,"Invalid vectors in opDot");
                assert(this.v2.dof.arrayStruct == this.v1.dof.arrayStruct,"Different array structs in opDot");
                assert(this.v1.dof.contiguous && this.v2.dof.contiguous,"not implemented for non contiguous arrays");
                if (this.v1.dof.length==0) return;
                auto d1=a2NA(this.v1.dof.support.data);
                auto d2=a2NA(this.v2.dof.support.data);
                auto rAtt=dot(d1,d2);
                addRes(cast(U)rAtt);
            } catch(Exception eTmp){
                this.e=new Exception("exception in doDofLoop",__FILE__,__LINE__,eTmp);
            }
        }
        void doOrientLoop(){
            try{
                if (this.e!is null) return;
                assert(this.v2.orient!is null,"Different vectors in opDot");
                assert(this.v2.orient.arrayStruct == this.v1.orient.arrayStruct,"Different array structs in opDot");
                assert(this.v1.orient.contiguous && this.v2.orient.contiguous,"not implemented for non contiguous arrays");
                if (this.v1.orient.length==0) return;
                auto d1=NArray!(T,2)([cast(index_type)4*T.sizeof,T.sizeof],[cast(index_type)this.v1.orient.support.length,3],
                    cast(index_type)0,(cast(T*)this.v1.orient.support.ptr)[0..4*this.v1.orient.support.length],0);
                auto d2=NArray!(T,2)([cast(index_type)4*T.sizeof,T.sizeof],[cast(index_type)this.v2.orient.support.length,3],
                    cast(index_type)0,(cast(T*)this.v2.orient.support.ptr)[0..4*this.v2.orient.support.length],0);
                auto rAtt=dotAll(d1,d2);
                addRes(cast(U)rAtt);
            } catch(Exception eTmp){
                this.e=new Exception("exception in doOrientLoop",__FILE__,__LINE__,eTmp);
            }
        }
        void doDot(){
            if (this.v1.pos!is null){
                Task("DynPVectorDotPos",&doPosLoop).autorelease.submitYield();
            } else {
                assert(this.v2.pos is null,"different vectors in dot");
            }
            if (this.v1.cell!is null){
                U resTmp=0;
                auto c1=this.v1.cell.h.cell[];
                auto c2=this.v2.cell.h.cell[];
                for(int i=0;i<9;++i){
                    resTmp+=c1[i]*c2[i];
                }
                addRes(resTmp);
            } else {
                assert(this.v2.cell is null,"different vectors in opDot");
            }
            if (this.v1.dof!is null){
                Task("DynPVectorDotDof",&doDofLoop).autorelease.submitYield();
            } else {
                assert(this.v2.dof is null,"different vectors in opDot");
            }
            if (this.v1.orient!is null){
                Task("DynPVectorDotOrient",&doOrientLoop).autorelease.submit();
            } else {
                assert(this.v2.orient is null,"different vectors in opDot");
            }
        }
    }
    U opDot(V,U=typeof(T.init+V.dtype.init))(V v2){ // could overlap the various loops...
        static assert(is(typeof(V.vGroup)),"V has to be a DynPVector");
        static assert(V.vGroup==group,"dot product only between same group");
        auto dOp=new DotOpStruct!(V,U); // keep on stack??? should be safe...
        dOp.v2=v2;
        dOp.v1=this;
        Task("DynPVectorDot",&dOp.doDot).autorelease.executeNow();
        if (dOp.e!is null) throw dOp.e;
        return dOp.res;
    }

    /// utility method for the squared euclidean (2-norm) of the vector: (this.dot(*this))
    T norm22(){
        return this.opDot(*this);
    }
    /// utility method for the euclidean (2-norm) of the vector: (this.dot(*this))
    T norm2(){
        return cast(T)sqrt(this.opDot(*this));
    }
    // this is *really* badly hacked in
    struct OpApplyStruct(int para){
        int res=0;
        Exception e=null;
        int delegate(ref T)loopBody;
        int delegate(ref size_t,ref T)loopBody2;
        DynPVector *pVect;
        
        void doPosLoop(){
            if (this.res!=0) return;
            try{
                static if (para!=0){
                    auto loopS=pVect.pos.pDataLoop;
                } else {
                    auto loopS=pVect.pos.sDataLoop;
                }
                auto resTmp=loopS.opApply(delegate int(ref Vector!(T,3) v){ // could be flattened using a NArray on pos.data like the dot...
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
        void doPosLoop2(){
            if (this.res!=0) return;
            try{
                static if (para!=0){
                    auto loopS=pVect.pos.pDataLoop;
                } else {
                    auto loopS=pVect.pos.sDataLoop;
                }
                auto resTmp=loopS.opApply(delegate int(ref size_t i,ref Vector!(T,3) v){ // could be flattened using a NArray on pos.data like the dot...
                    auto iB=3*i;
                    if (auto res=loopBody2(iB,v.x)) return res;
                    auto iB1=iB+1;
                    if (auto res=loopBody2(iB1,v.y)) return res;
                    auto iB2=iB+2;
                    if (auto res=loopBody2(iB2,v.z)) return res;
                    return 0;
                });
            } catch (Exception eTmp){
                e=new Exception("exception in doPosLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doDofLoop(){
            assert(pVect.dof.contiguous,"not implemented for non contiguous arrays");
            try{
                auto resTmp=pVect.dof.support.opApply(loopBody);
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doDofLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doDofLoop2(){
            assert(pVect.dof.contiguous,"not implemented for non contiguous arrays");
            try{
                size_t i0=pVect.pos.arrayStruct.dataLength *3;
                auto resTmp=pVect.dof.support.opApply(delegate int(ref size_t i,ref T el){
                    auto ii=i0+i;
                    return loopBody2(ii,el);
                });
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doDofLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doOrientLoop(){
            try{
                assert(pVect.orient.contiguous,"not implemented for non contiguous arrays");
                auto resTmp=pVect.orient.support.opApply(delegate int(ref Quaternion!(T) q){ // could be flattened using a NArray on pos.data...
                    if(auto r=loopBody(q.xyzw.x)) return r;
                    if(auto r=loopBody(q.xyzw.y)) return r;
                    if(auto r=loopBody(q.xyzw.z)) return r;
                    q.normalize();
                });
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doOrientLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doOrientLoop2(){
            try{
                size_t i0=3*pVect.pos.arrayStruct.dataLength+pVect.dof.arrayStruct.dataLength;
                assert(pVect.orient.contiguous,"not implemented for non contiguous arrays");
                auto resTmp=pVect.orient.support.opApply(delegate int(ref size_t i,ref Quaternion!(T) q){ // could be flattened using a NArray on pos.data...
                    auto iQ0=3*i+i0;
                    if(auto r=loopBody2(iQ0,q.xyzw.x)) return r;
                    auto iQ1=iQ0+1;
                    if(auto r=loopBody2(iQ1,q.xyzw.y)) return r;
                    auto iQ2=iQ0+2;
                    if(auto r=loopBody2(iQ2,q.xyzw.z)) return r;
                    q.normalize();
                });
                if (resTmp!=0) res=resTmp;
            } catch(Exception eTmp){
                e=new Exception("exception in doOrientLoop",__FILE__,__LINE__,eTmp);
                res=-1;
            }
        }
        void doOpApply(){
            if (pVect.pos!is null){
                static if (para!=0) Task("DynPVectorOpApplyPos",&doPosLoop).autorelease.submitYield();
                else doPosLoop();
            }
            if (pVect.dof!is null){
                static if (para!=0) Task("DynPVectorOpApplyDof",&doDofLoop).autorelease.submitYield();
                else doDofLoop();
            }
            if (pVect.orient!is null){
                static if (para!=0) Task("DynPVectorOpApplyOrient",&doOrientLoop).autorelease.submit();
                else doOrientLoop();
            }
            if (pVect.cell!is null){
                auto resTmp=pVect.cell.h.opApply(loopBody);
                if (resTmp!=0) res=resTmp;
            }
        }
        void doOpApply2(){
            if (pVect.pos!is null){
                static if (para!=0) Task("DynPVectorOpApplyPos2",&doPosLoop2).autorelease.submitYield();
                else doPosLoop2();
            }
            if (pVect.dof!is null){
                static if (para!=0) Task("DynPVectorOpApplyDof2",&doDofLoop2).autorelease.submitYield();
                else doDofLoop2();
            }
            if (pVect.orient!is null){
                static if (para!=0) Task("DynPVectorOpApplyOrient2",&doOrientLoop2).autorelease.submit();
                else doOrientLoop2();
            }
            if (pVect.cell!is null){
                auto i0=pVect.pos.arrayStruct.dataLength*3+pVect.orient.arrayStruct.dataLength*3
                    +pVect.dof.arrayStruct.dataLength;
                auto resTmp=pVect.cell.h.opApply(delegate int(ref size_t i, ref T el){
                    auto ii=i+i0;
                    return loopBody2(ii,el);
                });
                if (resTmp!=0) res=resTmp;
            }
        }
        int opApply(int delegate(ref T)loopBody){
            this.loopBody=loopBody;
            Task("DynPVectorOpApply",&this.doOpApply).autorelease.executeNow();
            if (this.e!is null) throw this.e;
            return this.res;
        }
        int opApply(int delegate(ref size_t,ref T)loopBody){
            this.loopBody2=loopBody;
            Task("DynPVectorOpApply",&this.doOpApply2).autorelease.executeNow();
            if (this.e!is null) throw this.e;
            return this.res;
        }
    }
    OpApplyStruct!(1) pLoop(size_t optimalBlockSize=1){
        OpApplyStruct!(1) oAp; // remove from stack???
        oAp.pVect=this;
        return oAp;
    }
    OpApplyStruct!(0) sLoop(){
        OpApplyStruct!(0) oAp; // remove from stack???
        oAp.pVect=this;
        return oAp;
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
    DynPVector!(V,group) dupT(V)(){
        DynPVector!(V,group) res;
        if (cell!is null) res.cell=cell.dupT!(V);
        if (pos!is null) res.pos=pos.dupT!(Vector!(V,3));
        if (orient!is null) res.orient=orient.dupT!(Quaternion!(V));
        if (dof!is null) res.dof=dof.dupT!(V);
        return res;
    }
    /// returns the convex hull of the kinds of all components
    KindRange gRange(){
        KindRange res;
        if (pos!is null){
            res=pos.kRange;
        }
        if (orient!is null){
            res=res.convexHull(orient.kRange);
        }
        if (dof!is null){
            res=res.convexHull(dof.kRange);
        }
        return res;
    }
    /// allocates an empty vector from the given pool
    static DynPVector allocFromPool(DynPVectStruct!(T) p,bool allocCell=false){
        DynPVector res;
        if (allocCell) res.cell=new Cell!(T)(Matrix!(T,3,3).identity,p.cellPeriod);
        res.pos=p.newPos();
        res.orient=p.newOrient();
        res.dof=p.newDof();
        return res;
    }
    mixin(serializeSome("dchem.sys.DynamicsVars("~T.stringof~")","cell|pos|orient|dof"));
    mixin printOut!();
}

/// matrix going from vectors in g1 to vectors in g2
class DynPMatrix(T,int g1,int g2){
    alias T dtype;
    enum{ rowGroupId=g1, colGroupId=g2 }
    NArray!(T,2) data;
    NArray!(T,2)[3][3] blocks;
    index_type[4] rowIdxs,colIdxs;
    DynPVectStruct!(T) rowGroup;
    DynPVectStruct!(T) colGroup;
    this(DynPVectStruct!(T)rowGroup,DynPVectStruct!(T)colGroup,NArray!(T,2) data=null){
        this.rowGroup=rowGroup;
        this.colGroup=colGroup;
        rowIdxs[0]=0;
        rowIdxs[1]=3*rowGroup.posStruct.dataLength;
        rowIdxs[2]=rowIdxs[1]+rowGroup.dofStruct.dataLength;
        rowIdxs[3]=rowIdxs[2]+4*rowGroup.orientStruct.dataLength;
        colIdxs[0]=0;
        colIdxs[1]=3*colGroup.posStruct.dataLength;
        colIdxs[2]=colIdxs[1]+colGroup.dofStruct.dataLength;
        colIdxs[3]=colIdxs[2]+4*colGroup.orientStruct.dataLength;
        if (data is null){
            this.data=empty!(T)([rowIdxs[3],colIdxs[3]],true); // fortran storage is more efficent for multiplication from the right
        } else {
            this.data=data;
            assert(data.shape[0]==rowIdxs[3],"invalid row size");
            assert(data.shape[1]==colIdxs[3],"invalid row size");
        }
        for (int iStripe=0;iStripe<3;++iStripe)
        for (int jStripe=0;jStripe<3;++jStripe){
            blocks[iStripe][jStripe]=data[Range(colIdxs[iStripe],colIdxs[iStripe+1]),Range(colIdxs[jStripe],colIdxs[jStripe+1])];
        }
    }
    
    /// matrix times the vector v1 in v2: v2 = scaleV2*v2+scaleRes*M*v1
    void matVectMult(V,int gg2,U,int gg1,U1=U,U2=U)(DynPVector!(V,gg2)v1,DynPVector!(U,gg1)v2,U1 scaleRes_=1,U2 scaleV2=0){
        static assert(gg1==g1 && gg2==g2,"mismatched groups");
        size_t[4] idxs;
        assert(v2.pos!is null && v2.orient!is null && v2.dof!is null,"target vector must be fully allocated");
        assert(rowIdxs[1]==3*v2.pos.length,"unexpected size of v2.pos");
        assert(rowIdxs[2]-rowIdxs[1]==v2.dof.length,"unexpected size of v2.dof");
        assert(rowIdxs[3]-rowIdxs[2]==4*v2.orient.length,"unexpected size of v2.orient");
        assert(0,"to do");
        assert(v2.pos.contiguous && v1.pos.contiguous,"not implemented for non contiguous segArrays");
        assert(v2.dof.contiguous && v1.dof.contiguous,"not implemented for non contiguous segArrays");
        assert(v2.orient.contiguous && v1.orient.contiguous,"not implemented for non contiguous segArrays");
        auto scaleRes=cast(U)scaleRes_;
        scope v2Pos=a2NA(v2.pos.support.basicData);
        scope v2Orient=a2NA(v2.orient.support.basicData);
        scope v2Dof=a2NA(v2.dof.support.basicData);
        U myscaleV2=cast(U)scaleV2;
        if (v1.pos!is null){
             assert(3*v1.pos.support.length==colIdxs[1],"unexpected size of v1.pos");
             scope v0=a2NA(v1.pos.support.basicData);
             dot(blocks[0][0],v0,v2Pos,scaleRes,myscaleV2);
             dot(blocks[1][0],v0,v2Dof,scaleRes,myscaleV2);
             dot(blocks[2][0],v0,v2Orient,scaleRes,myscaleV2);
             myscaleV2=1;
         }
         if (v1.dof!is null){
             assert(v1.dof.support.length==colIdxs[2]-colIdxs[1],"unexpected size of v1.pos");
             scope v0=a2NA(v1.dof.support.basicData);
             dot(blocks[0][2],v0,v2Pos,scaleRes,myscaleV2);
             dot(blocks[1][2],v0,v2Dof,scaleRes,myscaleV2);
             dot(blocks[2][2],v0,v2Orient,scaleRes,myscaleV2);
             myscaleV2=1;
         }
         if (v1.orient!is null){
             assert(3*v1.orient.support.length==colIdxs[3]-colIdxs[2],"unexpected size of v1.orient");
             scope v0=a2NA(v1.orient.support.basicData);
             dot(blocks[0][1],v0,v2Pos,scaleRes,myscaleV2);
             dot(blocks[1][1],v0,v2Dof,scaleRes,myscaleV2);
             dot(blocks[2][1],v0,v2Orient,scaleRes,myscaleV2);
             myscaleV2=1;
         }
        if (myscaleV2!=1){
            v2*=myscaleV2;
        }
    }
    /// transposed matrix times the vector v1 in v2: v2 = scaleV2*v2+scaleRes*M^T*v1
    void matTVectMult(V,int gg1,U,int gg2,U1=U,U2=U)(DynPVector!(V,gg1)v1,DynPVector!(U,gg2)v2,U1 scaleRes=1,U2 scaleV2=0){
        static assert(gg1==g1 && gg2==g2,"mismatched groups");
        size_t[4] idxs;
        assert(v2.pos!is null && v2.orient!is null && v2.dof!is null,"target vector must be fully allocated");
        assert(colIdxs[1]==3*v2.pos.length,"unexpected size of v2.pos");
        assert(colIdxs[2]-colIdxs[1]==v2.dof.length,"unexpected size of v2.dof");
        assert(colIdxs[3]-colIdxs[2]==4*v2.orient.length,"unexpected size of v2.orient");
        assert(v2.arrayMap.contiguous && v1.arrayMap.contiguous,"not implemented for non contiguous segArrays");
        scope v2Pos=a2NA(v2.pos.support.basicData);
        scope v2Orient=a2NA(v2.orient.support.basicData);
        scope v2Dof=a2NA(v2.dof.support.basicData);
        T myscaleV2=scaleV2;
        if (v1.pos!is null){
            assert(3*v1.pos.support.length==colIdxs[1],"unexpected size of v1.pos");
            scope v0=a2NA(v1.pos.support.basicData);
            dot(blocks[0][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[0][1],v0,v2Dof,scaleRes,myscaleV2,0,0);
            dot(blocks[0][2],v0,v2Orient,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (v1.dof!is null){
            assert(v1.dof.support.length==colIdxs[2]-colIdxs[1],"unexpected size of v1.dof");
            scope v0=a2NA(v1.dof.support.basicData);
            dot(blocks[1][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[1][1],v0,v2Dof,scaleRes,myscaleV2,0,0);
            dot(blocks[1][2],v0,v2Orient,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (v1.orient!is null){
            assert(3*v1.orient.support.length==colIdxs[3]-colIdxs[2],"unexpected size of v1.orient");
            scope v0=a2NA(v1.orient.support.basicData);
            dot(blocks[2][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[2][1],v0,v2Dof,scaleRes,myscaleV2,0,0);
            dot(blocks[2][2],v0,v2Orient,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (myscaleV2!=1){
            v2*=myscaleV2;
        }
    }
    /// outer multiplication
    void opOuter(V,int gg2,U,int gg1,U1=U,U2=U)(DynPVector!(V,gg2)v1,DynPVector!(U,gg1)v2,U1 scaleRes=1,U2 scaleThis=0){
        static assert(gg1==g1 && gg2==g2,"mismatched groups");
        assert(0,"to do"); // should be similar to dot, but using outer...
    }
    
    template dotOpRes(T,U,S...){
        static if(is(T==DynPMatrix)&& is(U.vGroup)){
            alias DynPVector!(T.dtype,T.colGroupId) dotOpRes;
        } else static if(is(U==DynPMatrix) && is(T.vGroup)){
            alias DynPVector!(T.dtype,T.rowGroupId) dotOpRes;
        } else {
            alias void dotOpRes;
        }
    }
    
    /// definition for the generic matrix vector multiplication
    static dotOpRes!(T,U,S) dotOp(T,U,S...)(T t,U u,S args){
        static if(is(T==DynPMatrix)&& is(U.vGroup)){
            static if(nArgs!(S)==0){
                auto res=t.colGroup.allocGroup!(g2)();
                t.matVectMult(u,res);
                return res;
            } else static if (is(typeof(args[0].vGroup))){
                t.matVectMult(u,args);
                return args[0];
            } else{
                auto res=t.colGroup.allocGroup!(g2)();
                t.matVectMult(u,res,args);
                return res;
            }
        } else static if(is(U==DynPMatrix) && is(T.vGroup)){
            static if(nArgs!(S)==0){
                auto res=u.colGroup.allocGroup!(g2)();
                u.matTVectMult(t,res);
                return res;
            } else static if (is(typeof(args[0].vGroup))){
                u.matVectMult(t,args);
                return args[0];
            } else{
                auto res=u.colGroup.allocGroup!(g2)();
                u.matVectMult(t,res,args);
                return res;
            }
        } else {
            static assert(0,"no dotOp for arguments:"~T.stringof~","~U.stringof~","~S.stringof);
        }
    }
    /// indexes, returns a copy...
    DynPVector!(T,g2)opIndex(size_t i){
        auto res=colGroup.allocGroup!(g2)();
        assert(0,"to do");
        return res;
    }
    /// definition for generic outer op
    static DynPMatrix outerOp(T,U,S...)(T v1,U v2,S args){
        static if (nArgs!(S)!=0 && is(S[0]==DynPMatrix)){
            args[0].opOuter(v1,v2,args[1..$]);
            return args[0];
        } else{
            //DynPMatrix res=new DynPMatrix(); // be smarter and take things from the vectors???
            //res.opOuter(v1,v2,args);
            assert(0,"to do");
            //return res;
        }
    }
    /// scalar product defined by this matrix: v1^T*M*v2
    /// could be optimized to avoid the use of a temporary
    T scalarProduct(V,U)(DynPVector!(V,g1)v1,DynPVector!(U,g2)v2){
        auto tmp=v1.emptyCopy;
        matVectMult(v2,tmp);
        assert(v1.arrayStruct==v2.arrayStruct,"redistribution of vectors not implemented");
        return v1.dot(v2.toGroup!(g1)());
    }
}

/// structure of the variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
class DynamicsVarsStruct(T){
    alias T dtype;

    /// structure (and pool) of the positions
    DynPVectStruct!(T) xGroup;
    /// structure (and pool) of the derivatives
    DynPVectStruct!(T) dxGroup;
    /// structure (and pool) of the dual derivatives (this is an extra separate group because if the overlap matrix
    /// is distributed the local vectors will indeed be incompatible)
    DynPVectStruct!(T) dualDxGroup;

    /// if the structures are compatible
    bool compatible(V)(DynamicsVarsStruct!(V) t){
        if (t is null) return false;
        return this is t || (xGroup.compatible(t.xGroup) && 
            dxGroup.compatible(t.dxGroup) && dualDxGroup.compatible(t.dualDxGroup));
    }
    override equals_t opEquals(Object o){
        auto t=cast(DynamicsVarsStruct)o;
        if (t is null) return 0;
        return this.compatible(t);
    }
    /// utility method, returns an empty position vector (like x)
    DynPVector!(T,XType) emptyX(bool allocCell=false){
        return DynPVector!(T,XType).allocFromPool(xGroup,allocCell);
    }
    /// utility method, returns an empty derivative vector (like dx,mddx)
    DynPVector!(T,DxType) emptyDx(bool allocCell=false){
        return DynPVector!(T,DxType).allocFromPool(dxGroup,allocCell);
    }
    /// utility method, returns an empty dual derivative vector
    DynPVector!(T,DualDxType) emptyDualDx(bool allocCell=false){
        return DynPVector!(T,DualDxType).allocFromPool(dualDxGroup,allocCell);
    }
    /// releases the cache of the groups
    void flush(){
        if (xGroup!is null) xGroup.flush;
        if (dxGroup!is null) dxGroup.flush;
        if (dualDxGroup!is null) dualDxGroup.flush;
    }
    /// releases the cache of the groups, and stops caching
    void stopCaching(){
        if (xGroup!is null) xGroup.stopCaching;
        if (dxGroup!is null) dxGroup.stopCaching;
        if (dualDxGroup!is null) dualDxGroup.stopCaching;
    }
    this(DynPVectStruct!(T)xGroup,DynPVectStruct!(T)dxGroup,DynPVectStruct!(T)dualDxGroup){
        this.xGroup=xGroup;
        this.dxGroup=dxGroup;
        this.dualDxGroup=dualDxGroup;
    }
    /// creates a copy that shares the pos/orient/dofStructs of another DynamicsVarsStruct
    static DynamicsVarsStruct copyFrom(V)(DynamicsVarsStruct!(V) s2){
        auto xGroup=new DynPVectStruct!(T)("xGroup",s2.xGroup.posStruct,s2.xGroup.orientStruct,s2.xGroup.dofStruct);
        auto dxGroup=new DynPVectStruct!(T)("dxGroup",s2.dxGroup.posStruct,s2.dxGroup.orientStruct,s2.dxGroup.dofStruct);
        auto dualDxGroup=new DynPVectStruct!(T)("dualDxGroup",s2.dualDxGroup.posStruct,s2.dualDxGroup.orientStruct,s2.dualDxGroup.dofStruct);
        auto res=new DynamicsVarsStruct(xGroup,dxGroup,dualDxGroup);
        static if(is(T==V)) {
            res.xGroup.allocPools(xGroup);
            res.dxGroup.allocPools(dxGroup);
            res.dualDxGroup.allocPools(dualDxGroup);
        } else {
            xGroup.allocPools();
            dxGroup.allocPools(xGroup);
            dualDxGroup.allocPools(dxGroup);
        }
        return res;
    }
    /// constructor, allocates the segmented array structures
    this(SubMapping fullSystem,KindRange kRange){
        xGroup=new DynPVectStruct!(T)("xGroup",fullSystem,kRange,SegmentedArrayStruct.Flags.None);
        dxGroup=new DynPVectStruct!(T)("dxGroup",fullSystem,kRange,SegmentedArrayStruct.Flags.None);
        dualDxGroup=new DynPVectStruct!(T)("dualDxGroup",fullSystem,kRange,SegmentedArrayStruct.Flags.None);
    }
    /// the structures have been defined and can be consolidated (if possible)
    /// sets up the pools
    void freezeStructs(){
        dxGroup.consolidateStructs(xGroup);
        dualDxGroup.consolidateStructs(dxGroup);
        xGroup.allocPools();
        dxGroup.allocPools(xGroup);
        dualDxGroup.allocPools(dxGroup);
    }
    
    /// checks that the given vector is allocated
    void checkAllocDynPVect(V,int i)(DynPVector!(V,i)* v){
        static if (i==0){
            auto pool=xGroup;
        } else static if (i==1){
            auto pool=dxGroup;
        } else static if (i==2){
            auto pool=dualDxGroup;
        } else {
            static assert(0,"unexpected vector group "~ctfe_i2a(i));
        }
        static if (is(V==T)){
            if (v.pos is null)    v.pos=pool.newPos();
            if (v.orient is null) v.orient=pool.newOrient();
            if (v.dof is null)    v.dof=pool.newDof();
        } else {
            if (v.pos is null)    {
                auto vMap=new SegArrMemMap!(Vector!(V,3))(pool.posStruct);
                v.pos=vMap.newArray();
            }
            if (v.orient is null) {
                auto oMap=new SegArrMemMap!(Vector!(V,3))(pool.orientStruct);
                v.orient=oMap.newArray();
            }
            if (v.dof is null)    {
                auto dofMap=new SegArrMemMap!(V)(pool.dofStruct);
                v.dof=dofMap.newArray();
            }
        }
    }
    
    /// returns a vector of the given type
    DynPVector!(V,i) dynPVector(int i)(){
        DynPVector!(V,i) res;
        checkAllocDynPVect(res);
        return res;
    }
}

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars(T){
    alias T dtype;

    /// structure of dynamic vars
    DynamicsVarsStruct!(T) dVarStruct;

    /// potential energy of the system (NAN if unknown)
    Real potentialEnergy;
    /// position vector
    DynPVector!(T,XType) x;
    /// velocities vector
    DynPVector!(T,DxType) dx;
    /// forces vector
    DynPVector!(T,DxType) mddx;
    
    bool isDummy(){
        return x.isDummy && dx.isDummy && mddx.isDummy;
    }
    /// returns a copy with a nullified cell
    DynamicsVars nullCell(){
        DynamicsVars res=*this;
        res.x.cell=null;
        res.dx.cell=null;
        res.mddx.cell=null;
        return res;
    }
    
    /// ensures that the positions are allocated
    void checkX(){
        dVarStruct.checkAllocDynPVect(&x);
    }
    /// ensures that the velocities are allocated
    void checkDx(){
        dVarStruct.checkAllocDynPVect(&dx);
    }
    /// ensures that the forces are allocated
    void checkMddx(){
        dVarStruct.checkAllocDynPVect(&mddx);
    }
    
    void opSliceAssignT(V)(ref DynamicsVars!(V) d2){
        assert(dVarStruct!is null);
        assert(dVarStruct.compatible(d2.dVarStruct),"incompatible structures");
        potentialEnergy=d2.potentialEnergy;
        x.opSliceAssignT!(V)(d2.x);
        dx.opSliceAssignT!(V)(d2.dx);
        mddx.opSliceAssignT!(V)(d2.mddx);
    }
    /// copies this dynvar to d2 (taking care of needed allocations)
    void copyTo(V)(ref DynamicsVars!(V) d2){
        if (!x.isDummy) {
            d2.checkX();
        } else if (!d2.x.isDummy) {
            d2.x.giveBack();
        }
        if (!dx.isDummy) {
            d2.checkDx();
        } else if (!d2.dx.isDummy) {
            d2.dx.giveBack();
        }
        if (!mddx.isDummy) {
            d2.checkMddx();
        } else if (!d2.mddx.isDummy) {
            ds.mddx.giveBack();
        }
        opSliceAssignT!(V)(d2);
    }
    void opSliceAssignEl(T val){
        assert((dVarStruct !is null),"invalid structures");
        potentialEnergy=0;
        x.opSliceAssignEl(val);
        dx.opSliceAssignEl(val);
        mddx.opSliceAssignEl(val);
    }
    void opSliceAssign(V)(V v){
        static if (is(typeof(opSliceAssignEl(v)))){
            opSliceAssignEl(v);
        } else static if (is(typeof(this.opSliceAssignT(v)))){
            this.opSliceAssignT(v);
        } else static if (is(typeof(v.copyTo(*this)))){
            v.copyTo(*this);
        } else {
            static assert(0,"cannot assign from "~V.stringof~" to DynamicsVars!("~T.stringof~")");
        }
    }
    void opBypax(V)(DynamicsVars!(V) v,V a,T b){
        x.opBypax(v.x,a,b);
        dx.opBypax(v.dx,a,b);
        mddx.opBypax(v.mddx,a,b);
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
    DynamicsVars!(V) dupT(V=T)(){
        DynamicsVars!(V) res;
        static if(is(T==V)){
            res.dVarStruct=dVarStruct;
        } else {
            res.dVarStruct=DynamicsVarsStruct!(V).copyFrom!(T)(this.dVarStruct);
        }
        res.potentialEnergy=potentialEnergy;
        res.x=x.dupT!(V);
        res.dx=dx.dupT!(V);
        res.mddx=mddx.dupT!(V);
        return res;
    }
    DynamicsVars dup(){
        return dupT!(T)();
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
    /// gives back the vector (memory might be reused)
    void giveBack(){
        x.giveBack();
        dx.giveBack();
        mddx.giveBack();
        clear();
    }
    /// clears the vector
    void clear(){
        potentialEnergy=T.init;
        x.clear;
        dx.clear;
        mddx.clear;
        dVarStruct=null;
    }
    /// deallocs the data stored in the vector
    void deallocData(){
        x.giveBack();
        dx.giveBack();
        mddx.giveBack();
        potentialEnergy=T.init;
    }
    mixin(serializeSome("dchem.sys.DynamicsVars("~T.stringof~")","potentialEnergy|x|dx|mddx"));
    mixin printOut!();
}
// to have minimal compiletime checks when compiling this module alone
version(InstantiateSome){
    pragma(msg,"DynVars InstantiateSome");
    private{
        DynPVectStruct!(double) a;
        DynPVector!(double,DxType) v;
        DynPMatrix!(double,DxType,DualDxType) m;
        DynamicsVars!(double) dv;
    }
}