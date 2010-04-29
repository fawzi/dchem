/// Particle system, the basic way to represent the state of a system
/// author: Fawzi
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
import blip.t.core.sync.Mutex;
import blip.narray.NArray;
import gobo.blas.Types:BlasTypeForType;
import blip.t.core.Traits: ctfe_i2a;

/// a pool for dynamical properties
/// setup as follow: init, update *Structs, possibly consolidateStructs, allocPools, possibly consolidate
class DynPPool(T){
    char[] name;
    SegmentedArrayStruct posStruct;
    SegmentedArrayStruct orientStruct;
    SegmentedArrayStruct dofStruct;
    
    SegArrPool!(Vector!(T,3)) poolPos;
    SegArrPool!(Quaternion!(T)) poolOrient;
    SegArrPool!(T) poolDof;
    
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
    void allocPools(DynPPool!(T) baseVal=null)
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
    void consolidate(DynPPool!(T) baseVal)
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
    void consolidateStructs(DynPPool!(T) baseVal)
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
}

/// matrix going from vectors in g1 to vectors in g2
class DynPMatrix(T,int g1,int g2){
    alias T dtype;
    enum{ rowGroupId=g1, colGroupId=g2 }
    NArray!(T,2) data;
    NArray!(T,2)[3][3] blocks;
    index_type[4] rowIdxs,colIdxs;
    DynPPool!(T) rowGroup;
    DynPPool!(T) colGroup;
    this(DynPPool!(T)rowGroup,DynPPool!(T)colGroup,NArray!(T,2) data=null){
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
}

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

/// simple interface that offers the basic derivative transfer types
/// this is nt meant to give a full differential geometry like setup, but to give enough
/// to implement the basic algorithms in a clean way without too much overheaad (thus not all the
/// guarantees of an exact treatement can be assumed)
/// one should make no assumptions about any regularity beyond continuity for large displacements
char[] DerivTransferTMixin(char[] t){
    return `
    /// adds to res in the tangential space of pos (derivative space) the equivalent (to first order)
    /// to the one that goes from pos.x to pos.x+diffP, or better pos.x-0.5*diffP to pos.x+0.5*diffP
    void addToTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,0)diffP,DynPVector!(`~t~`,1)res);
    /// adds to the vector in res the derivative deriv from the dual tangential space of pos.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order.
    void addFromDualTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,2)deriv,DynPVector!(`~t~`,0)res);
    /// transfers from a foreign dual tangential space, this is conceptually equivalent to parallel
    /// transfer, but needs to be valid only if x and pos2.x are close, in particular
    /// a direct transfer res+=derivAtPos2 is always allowed.
    /// uses the dual tangential space as it should be more "uniform"
    void addFromForeignDualTangentialSpace(ParticleSys!(`~t~`)pos1, ParticleSys!(`~t~`)pos2,DynPVector!(`~t~`,2)derivAtPos2,DynPVector!(`~t~`,2)res);
    `;
}
/// implements DerivTransferT using the templates that end with T (to work around compiler bugs in implementing interfaces with alias)
char[] derivTransferTMixin(char[]t){
    return `
    void addToTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,0)diffP,DynPVector!(`~t~`,1)res){
        addToTangentialSpaceT(pos,diffP,res);
    }
    void addFromDualTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,2)deriv,DynPVector!(`~t~`,0)res){
        addFromDualTangentialSpaceT!(`~t~`)(pos,deriv,res);
    }
    void addFromForeignDualTangentialSpace(ParticleSys!(`~t~`)pos1,ParticleSys!(`~t~`)pos2,DynPVector!(`~t~`,2)derivAtPos2,DynPVector!(`~t~`,2)res){
        addFromForeignDualTangentialSpaceT!(`~t~`)(pos1,pos2,derivAtPos2,res);
    }
    void setOverlap(DynPMatrix!(`~t~`,1,2)m){
        setOverlapT!(`~t~`)(m);
    }
    `;
}

/// interface that implements Real and LowP versions (multiple inheritane of templated interface seems to trigger bugs)
interface DerivTransfer{
    mixin(DerivTransferTMixin("Real")~DerivTransferTMixin("LowP"));
    /// should return true if the overlap is different from the identity
    bool specialOverlap();
    /// should set the overlap, matrix should use 2d distribution, at the moment it uses just pos,orient,dof sequence
    void setOverlap(DynPMatrix!(BlasTypeForType!(Real),1,2)m);
    static if(! is(BlasTypeForType!(Real)==BlasTypeForType!(LowP))){
        /// should set the overlap, matrix should use 2d distribution, at the moment it uses just pos,orient,dof sequence
        void setOverlap(DynPMatrix!(BlasTypeForType!(LowP),1,2)m);
    }
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
class ParticleKind: Serializable,CopiableObjectI,DerivTransfer{
    enum DerivMap{
        SimpleMap,  /// a simple map from position to derivatives is enough
        ComplexMap, /// a complex map has to be performed
    }
    char[] _name;
    char[] name(){ return _name; }
    void name(char[] nName){ _name=nName; }
    char[] _potential;
    char[] potential(){ return _potential; }
    void potential(char[] nPotential){ _potential=nPotential; }
    char[] symbol;
    
    /// describes the elements of a single particle in a segmented array
    /// in particular if the simple copy to build the derivative should be valid
    struct SegAElements{
        index_type nElements; /// number of elements in position array
        index_type nDelements; /// number of elements in derivative arrays
        index_type[] skipElements; /// ranges to skip in the array
        index_type[] skipDelements; /// ranges to skip in deriv array
        index_type[] el2DelMap; /// mapping of the contiguous sequences
        DerivMap derivMap=DerivMap.SimpleMap; /// how the derivatives are mapped
        Object _lockObj;
        mixin(serializeSome("dchem.SegAElements",`nElements: number of elements in position array
            nDelements: number of elements in derivative arrays
            skipElements: ranges to skip in the array
            skipDelements: ranges to skip derivatives array
            el2DelMap: mapping of the contiguous sequences
            derivMap: kind of mapping`));
        mixin printOut!();
        Object lockObj(){
            if (_lockObj !is null) return _lockObj;
            synchronized{
                if (_lockObj !is null){
                    _lockObj=new Mutex();
                }
            }
        }
        /// copies arr to darr
        void addArrayToDArray(T,V)(BulkArray!(T)arr,BulkArray!(V)darr,index_type blockSizeArr, index_type blockSizeDarr){
            assert(blockSizeArr>=nElements,"invalid block size in arr");
            assert(blockSizeDarr>=nDelements,"invalid block size in darr");
            if (derivMap==DerivMap.SimpleMap){
                assert(nElements==nDelements);
                alias typeof(arr.basicData[0]) ArrBData;
                scope a1=NArray!(ArrBData,2)([cast(index_type)(T.sizeof*blockSizeArr),ArrBData.sizeof],
                    [arr.basicData.length/(T.sizeof/ArrBData.sizeof),nElements],0,arr.basicData,0);
                alias typeof(darr.basicData[0]) DarrBData;
                scope a2=NArray!(DarrBData,2)([cast(index_type)(V.sizeof*blockSizeArr),DarrBData.sizeof],
                    [darr.basicData.length/(V.sizeof/DarrBData.sizeof),nDelements],0,darr.basicData,0);
                a2+=a1;
            } else {
                assert(el2DelMap.length>1&&(el2DelMap.length-2)%3==0,"unexpected el2DelMap length");
                index_type resEl=blockSizeArr-nElements;
                index_type resDel=blockSizeDarr-nDelements;
                size_t nSeg=el2DelMap.length/3;
                size_t nPart=arr.length/nElements;
                auto p=arr.ptr;
                auto dp=darr.ptr;
                for (size_t iPart=nPart;iPart!=0;--iPart){
                    p+=el2DelMap[0];
                    dp+=el2DelMap[1];
                    for(size_t iseg=0;iseg!=nSeg;++iseg){
                        for (index_type i=0;i!=el2DelMap[iseg*3+2];--i){
                            static if(is(typeof((*dp)+=*p))){
                                (*dp)+=(*p);
                            } else {
                                *dp=(*dp)+(*p);
                            }
                            ++dp;
                            ++p;
                        }
                        p+=el2DelMap[iseg*3+3];
                        dp+=el2DelMap[iseg*3+4];
                    }
                    p+=resEl;
                    dp+=resDel;
                }
            }
        }
        /// copies arr to darr
        void addDArrayToArray(T,V)(BulkArray!(V)darr,BulkArray!(T)arr,index_type blockSizeDarr,index_type blockSizeArr){
            assert(blockSizeArr>=nElements,"invalid block size in arr");
            assert(blockSizeDarr>=nDelements,"invalid block size in darr");
            if (derivMap==DerivMap.SimpleMap){
                assert(nElements==nDelements);
                alias typeof(darr.basicData[0]) DarrBData;
                scope a1=NArray!(DarrBData,2)([cast(index_type)(V.sizeof*blockSizeArr),DarrBData.sizeof],
                    [darr.basicData.length/(V.sizeof/DarrBData.sizeof),nDelements],0,darr.basicData,0);
                alias typeof(arr.basicData[0]) ArrBData;
                scope a2=NArray!(ArrBData,2)([cast(index_type)(T.sizeof*blockSizeArr),ArrBData.sizeof],
                    [arr.basicData.length/(T.sizeof/ArrBData.sizeof),nElements],0,arr.basicData,0);
                a2+=a1;
            } else {
                assert(el2DelMap.length>2&&(el2DelMap.length)%3==0,"unexpected el2DelMap length");
                index_type resEl=blockSizeArr-nElements;
                index_type resDel=blockSizeDarr-nDelements;
                size_t nSeg=el2DelMap.length/3;
                size_t nPart=arr.length/nElements;
                auto p=arr.ptr;
                auto dp=darr.ptr;
                for (size_t iPart=nPart;iPart!=0;--iPart){
                    for(size_t iseg=0;iseg!=nSeg;++iseg){
                        for (index_type i=0;i!=el2DelMap[iseg*3];--i){ // use the safer < instead of != ?
                            static if(is(typeof((*p)+=*dp))){
                                (*p)+= (*dp);
                            } else {
                                *p=(*p)+(*dp);
                            }
                            ++dp;
                            ++p;
                        }
                        p+=el2DelMap[iseg*3+1];
                        dp+=el2DelMap[iseg*3+2];
                    }
                    p+=resEl;
                    dp+=resDel;
                }
            }
        }
        /// sets overlap
        void setOverlap(T)(NArray!(T,2)overlap,index_type blockSizeArr,index_type blockSizeDarr){
            assert(blockSizeArr>=nElements,"invalid block size in arr");
            assert(blockSizeDarr>=nDelements,"invalid block size in darr");
            assert(overlap.shape[0]%blockSizeArr==0,"unexpected overlap size");
            assert(overlap.shape[1]%blockSizeDarr==0,"unexpected overlap size");
            if (derivMap==DerivMap.SimpleMap){
                diag(overlap)[]=cast(T)1;
            } else {
                assert(el2DelMap.length>2&&(el2DelMap.length)%3==0,"unexpected el2DelMap length");
                index_type resEl=blockSizeArr-nElements;
                index_type resDel=blockSizeDarr-nDelements;
                size_t nSeg=el2DelMap.length/3;
                index_type arrIdx=0;
                index_type darrIdx=0;
                while (arrIdx!=overlap.shape[1] && darrIdx!=overlap.shape[0]){
                    for(size_t iseg=0;iseg!=nSeg;++iseg){
                        for (index_type i=0;i!=el2DelMap[iseg*3];--i){ // use the safer < instead of != ?
                            overlap[arrIdx,darrIdx]=cast(T)1;
                            ++arrIdx;
                            ++darrIdx;
                        }
                        arrIdx+=el2DelMap[iseg*3+1];
                        darrIdx+=el2DelMap[iseg*3+2];
                    }
                    arrIdx+=resEl;
                    darrIdx+=resDel;
                }
                assert(arrIdx==overlap.shape[1] && darrIdx==overlap.shape[0]);
            }
        }
        /// rebuilds the el2DelMap and updates the derivMap
        void rebuildEl2DelMap(){
            synchronized(lockObj){
                if (skipElements.length==0 && skipDelements.length==0){
                    derivMap=DerivMap.SimpleMap;
                    el2DelMap=null;
                } else {
                    index_type elPos=0,skipElPos=0;
                    index_type delPos=0,skipDelPos=0;
                    index_type resPos=0;
                    while(elPos<nElements && delPos<nElements){
                        index_type nMax=nElements-elPos;
                        if(skipElPos<skipElements.length){
                            if (elPos<skipElements[skipElPos]){
                                nMax=skipElements[skipElPos]-elPos;
                            } else {
                                nMax=0;
                            }
                        }
                        if (skipDelPos<skipDelements.length){
                            if (delPos<skipDelements[skipDelPos]){
                                nMax=min(nMax,skipDelements[skipDelPos]);
                            } else {
                                nMax=0;
                            }
                        }
                        if (resPos>=el2DelMap.length){
                            el2DelMap~=[cast(index_type)0,0,0];
                        }
                        assert(el2DelMap.length>resPos+2,"internal error el2DelMap");
                        el2DelMap[resPos]=nMax;
                        elPos+=nMax;
                        ++resPos;
                        if(skipElPos<skipElements.length){
                            assert(skipElPos+1<skipElements.length);
                            if (elPos<=skipElements[skipElPos]){
                                assert(elPos==skipElements[skipElPos]);
                                el2DelMap[resPos]=skipElements[skipElPos+1]-elPos;
                                skipElPos+=2;
                            } else {
                                el2DelMap[resPos]=0;
                            }
                        } else {
                            el2DelMap[resPos]=0;
                        }
                        elPos+=el2DelMap[resPos];
                        ++resPos;
                        if (skipDelPos<skipDelements.length){
                            assert(skipDelPos+1<skipDelements.length);
                            if (delPos<=skipDelements[skipDelPos]){
                                assert(delPos==skipDelements[skipDelPos]);
                                el2DelMap[resPos]=skipDelements[skipDelPos+1]-delPos;
                                skipDelPos+=2;
                            } else {
                                el2DelMap[resPos]=0;
                            }
                        } else {
                            el2DelMap[resPos]=0;
                        }
                        delPos+=el2DelMap[resPos];
                        ++resPos;
                    }
                    assert(elPos==nElements && delPos==nDelements);
                }
            }
        }
        /// adds the given number of elements to the particles of this type.
        /// if complexMapping is true the simplified mapping between elements and derivatives is
        /// skipped (this is required if nEl!=nDel).
        /// returns the indexes of the first element and delement added to the particle with this
        /// operation
        void addElements(index_type nEl,index_type nDel,bool complexMapping ,
            out index_type elIdx, out index_type delIdx)
        {
            assert(nEl==nDel || complexMapping,"the simple mapping works only if nEl==nDel");
            synchronized(lockObj){
                elIdx=nElements;
                delIdx=nDelements;
                nElements+=nEl;
                nDelements+=nDel;
                if (complexMapping){
                    if (skipElements.length>0 && skipElements[$-1]==elIdx){
                        skipElements[$-1]=nElements;
                    } else {
                        skipElements~=[elIdx,nElements];
                    }
                    if (skipDelements.length>0 && skipDelements[$-1]==delIdx){
                        skipDelements[$-1]=nDelements;
                    } else {
                        skipDelements~=[delIdx,nDelements];
                    }
                }
                rebuildEl2DelMap();
            }
        }
    }
    size_t subParticles; /// number of subparticles for each particle of this type
    SegAElements posEls; /// per particle position like (Vector(T,3)) elements
    SegAElements orientEls; /// per particle orientation like (Quaternion(T)) elements
    SegAElements dofEls; /// per particle generic degrees of freedom (T) elements
    DerivTransfer[] transferHandlers; /// extra functions that are called to handle the complex part of the derivative transfers
    LevelIdx _level;
    LevelIdx level(){ return _level; }
    void level(LevelIdx l){ _level=l; }
    KindIdx pKind;
    
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(this))("ParticleKind");
        metaI.addFieldOfType!(char[])("name","name of this particle kind",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(ubyte)("level","level of the particle",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(ushort)("pKind","kind index");
        metaI.addFieldOfType!(SegAElements)("posEls","number of 3D space elements",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(SegAElements)("dofEls","number generalized degrees of freedom elements",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(SegAElements)("orientEls","number of oriantation elements",
            SerializationLevel.debugLevel);
        metaI.addFieldOfType!(size_t)("subParticles","number of subParticles",
            SerializationLevel.debugLevel);
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    /// just for internal use
    this(){}
    
    this(char[] pName,LevelIdx pLevel,KindIdx kindIdx,char[] symbol=null,char[] potential=null){
        this._name=pName;
        this.symbol=symbol;
        this._potential=potential;
        this._level=pLevel;
        this.pKind=kindIdx;
        posEls._lockObj=this;
        dofEls._lockObj=this;
        orientEls._lockObj=this;
        
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
        posEls=p.posEls;
        orientEls=p.orientEls;
        dofEls=p.dofEls;
        subParticles=p.subParticles;
    }
    
    void serial(S)(S s){
        s.field(metaI[0],_name);
        auto ui=cast(ubyte*)&_level;
        s.field(metaI[1],*ui);
        auto us=cast(ushort*)&pKind;
        s.field(metaI[2],*us);
        s.field(metaI[3],posEls);
        s.field(metaI[4],dofEls);
        s.field(metaI[5],orientEls);
        s.field(metaI[6],subParticles);
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
        foreach (pool;[&p.dynVars.xGroup]){
            if (pool.posStruct.addToKindDim(pKind,posEls.nElements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to posStruct for its kind for group "~pool.name);
            }
            if (pool.orientStruct.addToKindDim(pKind,orientEls.nElements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to orientStruct for its kind for group "~pool.name);
            }
            if (pool.dofStruct.addToKindDim(pKind,dofEls.nElements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to dofStruct for its kind for group "~pool.name);
            }
        }
        foreach (pool;[&p.dynVars.dxGroup,&p.dynVars.dualDxGroup]){
            if (pool.posStruct.addToKindDim(pKind,posEls.nDelements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to posStruct for its kind for group "~pool.name);
            }
            if (pool.orientStruct.addToKindDim(pKind,orientEls.nDelements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to orientStruct for its kind for group "~pool.name);
            }
            if (pool.dofStruct.addToKindDim(pKind,dofEls.nDelements)!=0){
                throw new Exception("particleKind '"~name~"' should be the first to add to dofStruct for its kind for group "~pool.name);
            }
        }
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
    
    /// adds to res in the tangential space of pos (derivative space) the equivalent (to first order)
    /// to the one that goes from pos.x-0.5*diffP to pos.x+0.5*diffP
    void addToTangentialSpaceT(T)(ParticleSys!(T)pos,DynPVector!(T,0)diffP,DynPVector!(T,1)res){
        if (pKind in res.pos.kRange && pKind in diffP.pos.kRange){
            posEls.addArrayToDArray(res.pos[pKind],res.pos[pKind],
                diffP.pos.arrayStruct.kindDim(pKind),res.pos.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.orient.kRange && pKind in diffP.orient.kRange){
            orientEls.addArrayToDArray(res.orient[pKind],res.orient[pKind],
                diffP.orient.arrayStruct.kindDim(pKind),res.orient.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.dof.kRange && pKind in diffP.dof.kRange){
            dofEls.addArrayToDArray(res.dof[pKind],res.dof[pKind],
                diffP.dof.arrayStruct.kindDim(pKind),res.dof.arrayStruct.kindDim(pKind));
        }
        foreach(obj;transferHandlers){
            obj.addToTangentialSpace(pos,diffP,res);
        }
    }
    /// adds to the vector in res the derivative deriv from the tangential space of pos.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order.
    void addFromDualTangentialSpaceT(T)(ParticleSys!(T)pos,DynPVector!(T,2)deriv,DynPVector!(T,0)res){
        if (pKind in res.pos.kRange && pKind in deriv.pos.kRange){
            posEls.addDArrayToArray(res.pos[pKind],res.pos[pKind],
                deriv.pos.arrayStruct.kindDim(pKind),res.pos.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.orient.kRange && pKind in deriv.orient.kRange){
            orientEls.addDArrayToArray(res.orient[pKind],res.orient[pKind],
                deriv.orient.arrayStruct.kindDim(pKind),res.orient.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.dof.kRange && pKind in deriv.dof.kRange){
            dofEls.addDArrayToArray(res.dof[pKind],res.dof[pKind],
                deriv.dof.arrayStruct.kindDim(pKind),res.dof.arrayStruct.kindDim(pKind));
        }
        foreach(obj;transferHandlers){
            obj.addFromDualTangentialSpace(pos,deriv,res);
        }
    }
    /// transfers from a foreign tangential space, this is conceptually equivalent to parallel
    /// transfer, but needs to be valid only if x and pos2.x are close, in particular
    /// a direct transfer res+=derivAtPos2 is always allowed (and is what is done by default)
    void addFromForeignDualTangentialSpaceT(T)(ParticleSys!(T)pos1,ParticleSys!(T)pos2,
        DynPVector!(T,2)derivAtPos2,DynPVector!(T,2)res){
        res.axpby(derivAtPos2,1,1);
        foreach(obj;transferHandlers){
            obj.addFromForeignDualTangentialSpace(pos1,pos2,derivAtPos2,res);
        }
    }
    /// should return true if the overlap is different from the identity
    bool specialOverlap(){
        foreach(obj;transferHandlers){
            if (obj.specialOverlap) return true;
        }
        return false;
    }
    /// should set the overlap, matrix should use 2d distribution
    void setOverlapT(T)(DynPMatrix!(T,1,2)m){
        assert(m.rowGroup.posStruct==m.colGroup.posStruct && m.rowGroup.orientStruct==m.colGroup.orientStruct &&
            m.rowGroup.dofStruct==m.colGroup.dofStruct,"only col=row implemented");
        if (pKind in m.rowGroup.gRange && pKind in m.colGroup.gRange){
            if (pKind in m.rowGroup.posStruct.kRange && pKind in m.colGroup.posStruct.kRange){
                auto rKStart=pKind-m.rowGroup.posStruct.kRange.kStart;
                auto rKStarts=m.rowGroup.posStruct.kindStarts;
                auto cKStart=pKind-m.colGroup.posStruct.kRange.kStart;
                auto cKStarts=m.colGroup.posStruct.kindStarts;
                posEls.setOverlap(m.blocks[0][0][Range(rKStarts[rKStart],rKStarts[rKStart+1]),
                    Range(cKStarts[cKStart],cKStarts[cKStart+1])],m.rowGroup.posStruct.kindDim(pKind),
                    m.colGroup.posStruct.kindDim(pKind));
                rKStart=pKind-m.rowGroup.orientStruct.kRange.kStart;
                rKStarts=m.rowGroup.orientStruct.kindStarts;
                cKStart=pKind-m.colGroup.orientStruct.kRange.kStart;
                cKStarts=m.colGroup.orientStruct.kindStarts;
                orientEls.setOverlap(m.blocks[1][1][Range(rKStarts[rKStart],rKStarts[rKStart+1]),
                    Range(cKStarts[cKStart],cKStarts[cKStart+1])],m.rowGroup.orientStruct.kindDim(pKind),
                    m.colGroup.orientStruct.kindDim(pKind));
                rKStart=pKind-m.rowGroup.dofStruct.kRange.kStart;
                rKStarts=m.rowGroup.dofStruct.kindStarts;
                cKStart=pKind-m.colGroup.dofStruct.kRange.kStart;
                cKStarts=m.colGroup.dofStruct.kindStarts;
                dofEls.setOverlap(m.blocks[1][1][Range(rKStarts[rKStart],rKStarts[rKStart+1]),
                    Range(cKStarts[cKStart],cKStarts[cKStart+1])],m.rowGroup.dofStruct.kindDim(pKind),
                    m.colGroup.dofStruct.kindDim(pKind));
            }
        }
        foreach(obj;transferHandlers){
            obj.setOverlap(m);
        }
    }
    mixin(derivTransferTMixin("Real"));
    mixin(derivTransferTMixin("LowP"));
}

/// kind of the particle that repesents the whole system
class SysKind: ParticleKind{
    override char[] name(){ return "_SYSTEM_"; }
    override char[] potential() { return ""; }
    this(){}
    
    this(LevelIdx pLevel,KindIdx kindIdx){
        super("_SYSTEM_",pLevel,kindIdx,"","");
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
    
    U opDot(V,U=typeof(T.init+V.init))(DynPVector!(V,group) v2){ // could overlap the various loops...
        U res=0;
        Exception e=null;
        void doPosLoop(){
            try{
                if (e!is null) return;
                assert(v2.pos!is null,"Different vectors in opDot");
                assert(v2.pos.arrayStruct is pos.arrayStruct,"Different array structs in opDot");
                auto d1=a2NA(pos.data.basicData);
                auto d2=a2NA(v2.pos.data.basicData);
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
                auto d1=a2NA(dof.data.data);
                auto d2=a2NA(v2.dof.data);
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
    static DynPVector allocFromPool(DynPPool!(T) p,bool allocCell=false){
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

/// variables that normally an integrator should care about
// (put here mainly for tidiness reasons)
struct DynamicsVars(T){
    alias T dtype;
    Real potentialEnergy; /// potential energy of the system (NAN if unknown)

    /// structure (and pool) of the positions
    DynPPool!(T) xGroup;
    /// structure (and pool) of the derivatives
    DynPPool!(T) dxGroup;
    /// structure (and pool) of the dual derivatives (this is an extra separate group because if the overlap matrix
    /// is distributed the local vectors will indeed be incompatible)
    DynPPool!(T) dualDxGroup;
    
    /// position vector
    DynPVector!(T,0) x;
    /// velocities vector
    DynPVector!(T,1) dx;
    /// forces vector
    DynPVector!(T,1) mddx;
    
    /// reallocates the segmented array structures, invalidates everything
    void reallocStructs(SysStruct sys){
        if (xGroup!is null) xGroup.flush; // call stopCaching ?
        if (dxGroup!is null) dxGroup.flush; // call stopCaching ?
        if (dualDxGroup!is null) dualDxGroup.flush; // call stopCaching ?
        xGroup=new DynPPool!(T)("xGroup",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            SegmentedArrayStruct.Flags.None);
        dxGroup=new DynPPool!(T)("dxGroup",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            SegmentedArrayStruct.Flags.None);
        dualDxGroup=new DynPPool!(T)("dualDxGroup",sys.fullSystem,KindRange(sys.levels[0].kStart,sys.levels[$-1].kEnd),
            SegmentedArrayStruct.Flags.None);
        x.clear();
        dx.clear();
        mddx.clear();
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
    alias BlasTypeForType!(T) dtypeBlas;
    char[] name; /// name of the particle system
    ulong iteration;
    NotificationCenter nCenter;
    SysStruct sysStruct;
    
    DynamicsVars!(T) dynVars; /// dynamics variables
    HiddenVars hVars; /// possibly some hidden degrees of freedom
    struct DerivOverlap{
        DynPMatrix!(dtypeBlas,1,2) overlap; /// overlap of the derivatives (metric), might be null or outdated, access through maybeOverlap
        DynPMatrix!(dtypeBlas,2,1) overlapInv; /// pseudo inverse of the overlap, might be null or outdated, access through maybeOverlapInv
        bool specialOverlap=true;
        bool overlapUpToDate=false;
        bool overlapInvUpToDate=false;
        /// full reset of the structure
        void clear(){
            overlap=null;
            overlapInv=null;
            specialOverlap=true;
            overlapUpToDate=false;
            overlapInvUpToDate=false;
        }
        /// signals that the overlap has become invalid
        void invalidate(){
            overlapUpToDate=false;
            overlapInvUpToDate=false;
        }
        /// update of the overlap if specialOverlap && !overlapUpToDate
        void updateOverlap(ParticleSys pSys){
            synchronized(pSys){
                if (specialOverlap){
                    specialOverlap=pSys.specialDerivOverlap();
                }
                if (specialOverlap){
                    if (overlap is null || overlap.rowGroup.posStruct !is pSys.dynVars.dualDxGroup.posStruct
                        || overlap.rowGroup.orientStruct !is pSys.dynVars.dualDxGroup.orientStruct
                        || overlap.rowGroup.dofStruct !is pSys.dynVars.dualDxGroup.dofStruct
                        || overlap.colGroup.posStruct !is pSys.dynVars.dxGroup.posStruct
                        || overlap.colGroup.orientStruct !is pSys.dynVars.dxGroup.orientStruct
                        || overlap.colGroup.dofStruct !is pSys.dynVars.dxGroup.dofStruct)
                    {
                        auto rowGroup=pSys.dynVars.dxGroup;
                        auto colGroup=pSys.dynVars.dualDxGroup;
                        static if (! is(dtype==dtypeBlas)){
                            auto rGroup=new DynPPool!(dtypeBlas)("blasDxGroup",rowGroup.posStruct,
                                rowGroup.orientStruct,rowGroup.dofStruct);
                            auto cGroup=new DynPPool!(dtypeBlas)("blasDualDxGroup",colGroup.posStruct,
                                colGroup.orientStruct,colGroup.dofStruct);
                            overlap=new DynPMatrix!(dtypeBlas,1,2)(rGroup,cGroup);
                        } else {
                            overlap=new DynPMatrix!(dtypeBlas,1,2)(rowGroup,colGroup);
                        }
                    } else {
                        overlap.data[]=cast(dtypeBlas)0;
                    }
                    pSys.setDerivOverlap(overlap);
                }
                overlapUpToDate=true;
            }
        }
        /// update the inverse of the overlap
        void updateOverlapInv(ParticleSys pSys){
            synchronized(pSys){
                updateOverlap(pSys);
                if (specialOverlap){
                    if (overlapInvUpToDate) return;
                    assert(overlap!is null);
                    if (overlapInv is null || overlapInv.rowGroup!is pSys.dynVars.dxGroup 
                        || overlapInv.colGroup!is pSys.dynVars.dualDxGroup){
                        overlapInv=new DynPMatrix!(dtypeBlas,2,1)(overlap.colGroup,overlap.rowGroup);
                    }
                    // well inversion is a global thing, so most likely n!=m makes no sense, anyway...
                    index_type m=overlapInv.data.shape[0],n=overlapInv.data.shape[1],mn=min(n,m);
                    scope u=empty!(dtypeBlas)([m,mn]);
                    scope vt=empty!(dtypeBlas)([mn,n]);
                    scope s=empty!(RealTypeOf!(dtypeBlas))(mn);
                    scope s2=svd(overlap.data,u,s,vt); /// should use the 'O' method and drop vt
                    unaryOpStr!(`
                        if (abs(*aPtr0)>1.e-10) {
                            *aPtr0=1/(*aPtr0);
                        }
                    `,1,dtypeBlas)(s2);
                    if (n<=m){
                        vt*=repeat(s2,n,-1);
                        dot(u[Range(0,-1),Range(0,n)],vt,overlapInv.data);
                    } else {
                        u*=repeat(s2,m,0);
                        dot(u,vt[Range(0,m)],overlapInv.data);
                    }
                    overlapInvUpToDate=true;
                }
            }
        }
    }
    DerivOverlap derivOverlap;
    
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
    /// reallocates array structs
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
    
    /// kinds (actually the whole sys struct) have been created and are valid
    /// (can be used to add stuff to particles: extra dofs,...)
    void pKindsInitialSetup(){
        if (nCenter!is null)
            nCenter.notify("kindsSet",Variant(this));
        dynVars.freezeStructs();
    }
    /// system structure changed (particle added/removed, kinds added/removed)
    /// the segmented array structs should be initialized, and modifiable.
    /// positions,... are not yet available
    void sysStructChanged(){
        derivOverlap.clear();
        this.reallocStructs();
        foreach(pKind;sysStruct.particleKinds.data){
            pKind.sysStructChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("sysStructChanged",Variant(this));
        dynVars.freezeStructs();
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(){
        derivOverlap.invalidate();
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
    
    /// returns a vector in the tangential space (derivative space) equivalent (to first order)
    /// to the one that goes from x-0.5*scale*diffP to x+0.5*scale*diffP
    void addToTSpace(T)(DynPVector!(T,0)diffP,DynPVector!(T,1)res){
        foreach(k;diffP.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addToTangentialSpace(this,diffP,res);
        }
    }
    /// adds to the vector in res the derivative deriv.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order in particular
    /// r[]=0; r2[]=0; addFromTangentialSpace(deriv,r); r2[]=r-x; addToTangentialSpace(r2,d);
    /// implies that d and deriv might differ, but should become equal for small derivatives
    /// (in first order).
    void addFromDualTSpace(T)(DynPVector!(T,2)deriv,DynPVector!(T,0)res){
        foreach(k;deriv.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addFromTangentialSpace(this,deriv,res);
        }
    }
    /// transfers from a foreign tangential space, this is conceptually equivalent to parallel
    /// transfer, but needs to be valid only if x and pos2.x are close, in particular
    /// a direct transfer res+=derivAtPos2 is always allowed.
    void addFromForeignDualTSpace(T)(DynPVector!(T,0)pos2,DynPVector!(T,2)derivAtPos2,DynPVector!(T,2)res){
        foreach(k;derivAtPos2.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addFromForeignDualTangentialSpace(this,pos2,derivAtPos2,res);
        }
    }
    /// returns either an up to date overlap of the derivatives or null (which means that it is the identity)
    DynPMatrix!(dtypeBlas,1,2) maybeDerivOverlap(){
        if (!derivOverlap.specialOverlap) return null;
        if (derivOverlap.overlapUpToDate) return derivOverlap.overlap;
        derivOverlap.updateOverlap(this);
        if (!derivOverlap.specialOverlap) return null;
        return derivOverlap.overlap;
    }
    /// returns either an up to date inverse of the derivatives overlap or null (which means that it is the identity)
    DynPMatrix!(dtypeBlas,2,1) maybeDerivOverlapInv(){
        if (!derivOverlap.specialOverlap) return null;
        if (derivOverlap.overlapInvUpToDate) return derivOverlap.overlapInv;
        derivOverlap.updateOverlapInv(this);
        if (!derivOverlap.specialOverlap) return null;
        return derivOverlap.overlapInv;
    }
    bool specialDerivOverlap(){
        bool specOv=false;
        foreach(pKind;sysStruct.particleKinds.data){
            if (pKind.specialOverlap()){
                return true;
            }
        }
        return false;
    }
    void setDerivOverlap(DynPMatrix!(dtypeBlas,1,2) overlap){
        foreach(pKind;sysStruct.particleKinds.data){
            pKind.setOverlap(overlap);
        }
    }
    
    /// goes to the dual tangential space (multiplication with S^-1) and adds the result to dualDeriv
    void toDualTSpace(DynPVector!(T,1)deriv,DynPVector!(T,2)dualDeriv,T scaleRes=1,T scaleDualDeriv=0){
        auto overlapInv=maybeDerivOverlapInv();
        if (overlapInv!is null){
            overlapInv.matVectMult!(T,T)(deriv,dualDeriv,scaleRes,scaleDualDeriv);
        } else { // assumes two groups ar actually equal
            dualDeriv.axpby(deriv.toGroup!(2)(),scaleRes,scaleDualDeriv);
        }
    }
    /// comes back from the dual tangential space (multiplication with S)
    void fromDualTSpace(DynPVector!(T,2)dualDeriv,DynPVector!(T,1)deriv,T scaleRes=1,T scaleDeriv=0){
        auto overlap=maybeDerivOverlap();
        if (overlap!is null){
            overlap.matVectMult!(T,T)(dualDeriv,deriv,scaleRes,scaleDeriv);
        } else {
            deriv.axpby(dualDeriv.toGroup!(1)(),scaleRes,scaleDeriv);
        }
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
