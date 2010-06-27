/// dynamic variables of a system of particles, vectors, matrixes, pool, and a structure to keep all of them togheter
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

enum{
    XType=0,
    DxType=1,
    DualDxType=2,
}
/// this is at the same time a pool and a structure for a DynPVector
/// setup as follow: init, update *Structs, possibly consolidateStructs, allocPools, possibly consolidate
class DynPVectStruct(T){
    char[] name;
    int[3] cellPeriod;
    SegmentedArrayStruct posStruct;
    SegmentedArrayStruct orientStruct;
    SegmentedArrayStruct dofStruct;
    
    SegArrPool!(Vector!(T,3)) poolPos;
    SegArrPool!(Quaternion!(T)) poolOrient;
    SegArrPool!(T) poolDof;
    
    /// equal if it has the same structures
    override equals_t opEquals(Object o){
        auto t=cast(DynPVectStruct)o;
        if (t is null) return 0;
        return (posStruct==t.posStruct && orientStruct==t.orientStruct && dofStruct==t.dofStruct);
    }
    
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
            poolPos=new SegArrPool!(Vector!(T,3))(posStruct);
        }
        if (baseVal!is null && baseVal.poolOrient!is null && baseVal.poolOrient.arrayStruct is orientStruct){
            poolOrient=baseVal.poolOrient;
        } else {
            poolOrient=new SegArrPool!(Quaternion!(T))(orientStruct);
        }
        if (baseVal!is null && baseVal.poolDof!is null && baseVal.poolDof.arrayStruct is dofStruct){
            poolDof=baseVal.poolDof;
        } else {
            poolDof=new SegArrPool!(T)(dofStruct);
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

/// perform an operation on the segmented arrays and cell of a dynPVector
/// assumes the existence of a boolean variable named "weak" that if true suppress the exception
/// when some of the arrays are null and other aren't
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
            if ((!weak) && (`;
            foreach (i,n; namesLocal){
                if (i!=0) res~="||";
                res~=n;
                res~=pp;
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
struct DynPVector(T,int group){
    alias T dtype;
    enum{ vGroup=group }

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
    void axpby(V)(V x,T a=1,T b=1){
        static assert(is(V==DynPVector!(V.dtype,group)),"axpby only between DynPVectors of the same group, not "~V.stringof);
        enum{ weak=false }
        if ((cell is x.cell) && cell !is null){
            throw new Exception("identical cells in axpby",__FILE__,__LINE__);
        }
        auto y=this;
        mixin(dynPVectorOp(["x","y"],"y.axpby(x,a,b);",true,false));
    }
    void opMulAssign()(T scale){
        enum{ weak=false }
        auto x=this;
        mixin(dynPVectorOp(["x"],"x*=scale;",true,false));
    }
    void opMulAssign(V)(DynPVector!(V,group) y){
        enum{ weak=false }
        auto x=this;
        mixin(dynPVectorOp(["x","y"],"x*=y;",true,false));
    }
    
    void opSliceAssign(V)(DynPVector!(V,group) b){
        enum{ weak=false }
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
    
    U opDot(V,U=typeof(T.init+V.dtype.init))(V v2){ // could overlap the various loops...
        static assert(is(typeof(V.vGroup)),"V has to be a DynPVector");
        static assert(V.vGroup==group,"dot product only between same group");
        U res=0;
        Exception e=null;
        void addRes(U inc){
            static if(is(typeof(atomicAdd(res,inc)))){
                atomicAdd(res,inc);
            } else {
                synchronized{
                    res+=inc;
                }
            }
        }
        void doPosLoop(){
            try{
                if (e!is null) return;
                assert(v2.pos!is null,"Different vectors in opDot");
                assert(v2.pos.arrayStruct is pos.arrayStruct,"Different array structs in opDot");
                auto d1=a2NA(pos.data.basicData);
                auto d2=a2NA(v2.pos.data.basicData);
                auto rAtt=dot(d1,d2);
                addRes(cast(U)rAtt);
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
                auto d1=a2NA(dof.data.data);
                auto d2=a2NA(v2.dof.data.data);
                auto rAtt=dot(d1,d2);
                addRes(cast(U)rAtt);
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
                    cast(index_type)0,(cast(T*)orient.data.ptr)[0..4*orient.data.length],0);
                auto d2=NArray!(T,2)([cast(index_type)4*T.sizeof,T.sizeof],[cast(index_type)v2.orient.data.length,3],
                    cast(index_type)0,(cast(T*)v2.orient.data.ptr)[0..4*v2.orient.data.length],0);
                auto rAtt=dotAll(d1,d2);
                addRes(cast(U)rAtt);
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
                    addRes(resTmp);
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
        return res;
    }

    /// utility method for the squared euclidean (2-norm) of the vector: (this.dot(*this))
    T norm22(){
        return this.opDot(*this);
    }
    /// utility method for the euclidean (2-norm) of the vector: (this.dot(*this))
    T norm2(){
        return sqrt(this.opDot(*this));
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
        if (allocCell) res.cell=new Cell!(T)(Matrix!(T,3,3).identity,[0,0,0]);
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
        rowIdxs[2]=rowIdxs[1]+4*rowGroup.orientStruct.dataLength;
        rowIdxs[3]=rowIdxs[2]+rowGroup.dofStruct.dataLength;
        colIdxs[0]=0;
        colIdxs[1]=3*colGroup.posStruct.dataLength;
        colIdxs[2]=colIdxs[1]+4*colGroup.orientStruct.dataLength;
        colIdxs[3]=colIdxs[2]+colGroup.dofStruct.dataLength;
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
    void matVectMult(V,U)(DynPVector!(V,g2)v1,DynPVector!(U,g1)v2,U scaleRes=1,U scaleV2=0){
        size_t[4] idxs;
        assert(v2.pos!is null && v2.orient!is null && v2.dof!is null,"target vector must be fully allocated");
        assert(rowIdxs[1]==3*v2.pos.length,"unexpected size of v2.pos");
        assert(rowIdxs[2]-rowIdxs[1]==4*v2.orient.length,"unexpected size of v2.orient");
        assert(rowIdxs[3]-rowIdxs[2]==v2.dof.length,"unexpected size of v2.dof");
        scope v2Pos=a2NA(v2.pos.data.basicData);
        scope v2Orient=a2NA(v2.orient.data.basicData);
        scope v2Dof=a2NA(v2.dof.data.basicData);
        T myscaleV2=scaleV2;
        if (v1.pos!is null){
             assert(3*v1.pos.data.length==colIdxs[1],"unexpected size of v1.pos");
             scope v0=a2NA(v1.pos.data.basicData);
             dot!(T,2,V,1,U,1)(blocks[0][0],v0,v2Pos,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[1][0],v0,v2Orient,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[2][0],v0,v2Dof,scaleRes,myscaleV2);
             myscaleV2=1;
         }
         if (v1.orient!is null){
             assert(3*v1.orient.data.length==colIdxs[2]-colIdxs[1],"unexpected size of v1.orient");
             scope v0=a2NA(v1.orient.data.basicData);
             dot!(T,2,V,1,U,1)(blocks[0][1],v0,v2Pos,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[1][1],v0,v2Orient,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[2][1],v0,v2Dof,scaleRes,myscaleV2);
             myscaleV2=1;
         }
         if (v1.pos!is null){
             assert(v1.dof.data.length==colIdxs[3]-colIdxs[2],"unexpected size of v1.pos");
             scope v0=a2NA(v1.dof.data.basicData);
             dot!(T,2,V,1,U,1)(blocks[0][2],v0,v2Pos,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[1][2],v0,v2Orient,scaleRes,myscaleV2);
             dot!(T,2,V,1,U,1)(blocks[2][2],v0,v2Dof,scaleRes,myscaleV2);
             myscaleV2=1;
         }
        if (myscaleV2!=1){
            v2*=myscaleV2;
        }
    }
    /// transposed matrix times the vector v1 in v2: v2 = scaleV2*v2+scaleRes*M^T*v1
    void matTVectMult(V,U)(DynPVector!(V,g1)v1,DynPVector!(U,g2)v2,U scaleRes=1,U scaleV2=0){
        size_t[4] idxs;
        assert(v2.pos!is null && v2.orient!is null && v2.dof!is null,"target vector must be fully allocated");
        assert(colIdxs[1]==3*v2.pos.length,"unexpected size of v2.pos");
        assert(colIdxs[2]-colIdxs[1]==4*v2.orient.length,"unexpected size of v2.orient");
        assert(colIdxs[3]-colIdxs[2]==v2.dof.length,"unexpected size of v2.dof");
        scope v2Pos=a2NA(v2.pos.data.basicData);
        scope v2Orient=a2NA(v2.orient.data.basicData);
        scope v2Dof=a2NA(v2.dof.data.basicData);
        T myscaleV2=scaleV2;
        if (v1.pos!is null){
            assert(3*v1.pos.data.length==colIdxs[1],"unexpected size of v1.pos");
            scope v0=a2NA(v1.pos.data.basicData);
            dot(blocks[0][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[0][1],v0,v2Orient,scaleRes,myscaleV2,0,0);
            dot(blocks[0][2],v0,v2Dof,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (v1.orient!is null){
            assert(3*v1.orient.data.length==colIdxs[2]-colIdxs[1],"unexpected size of v1.orient");
            scope v0=a2NA(v1.orient.data.basicData);
            dot(blocks[1][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[1][1],v0,v2Orient,scaleRes,myscaleV2,0,0);
            dot(blocks[1][2],v0,v2Dof,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (v1.pos!is null){
            assert(v1.dof.data.length==colIdxs[3]-colIdxs[2],"unexpected size of v1.pos");
            scope v0=a2NA(v1.dof.data.basicData);
            dot(blocks[2][0],v0,v2Pos,scaleRes,myscaleV2,0,0);
            dot(blocks[2][1],v0,v2Orient,scaleRes,myscaleV2,0,0);
            dot(blocks[2][2],v0,v2Dof,scaleRes,myscaleV2,0,0);
            myscaleV2=1;
        }
        if (myscaleV2!=1){
            v2*=myscaleV2;
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

    override equals_t opEquals(Object o){
        auto t=cast(DynamicsVarsStruct)o;
        if (t is null) return 0;
        return this is t || (xGroup==t.xGroup && dxGroup==t.dxGroup && dualDxGroup==t.dualDxGroup);
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
        return new DynamicsVarsStruct(xGroup,dxGroup,dualDxGroup);
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
            if (v.pos is null)    v.pos=new SegmentedArray!(Vector!(V,3))(pool.posStruct);
            if (v.orient is null) v.orient=new SegmentedArray!(Quaternion!(V))(pool.orientStruct);
            if (v.dof is null)    v.dof=new SegmentedArray!(V)(pool.dofStruct);
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
    
    void opSliceAssign(V)(ref DynamicsVars!(V) d2){
        assert((dVarStruct is null && d2.dVarStruct is null)|| dVarStruct is d2.dVarStruct || dVarStruct == d2.dVarStruct,
            "incompatible structures");
        potentialEnergy=d2.potentialEnergy;
        x[]=d2.x;
        dx[]=d2.dx;
        mddx[]=d2.mddx;
    }
    void opSliceAssign()(T val){
        assert((dVarStruct is null && val.dVarStruct is null)|| dVarStruct is val.dVarStruct || dVarStruct == val.dVarStruct,
            "incompatible structures");
        potentialEnergy=0;
        x[]=val;
        dx[]=val;
        mddx[]=val;
    }
    void axpby(V)(DynamicsVars!(V) v,V a,T b){
        x.axpby(v.x,a,b);
        dx.axpby(v.dx,a,b);
        mddx.axpby(v.mddx,a,b);
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