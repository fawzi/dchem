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
import blip.core.Variant;
import blip.container.BitArray;
import blip.container.Deque;
import blip.container.BulkArray;
import blip.parallel.smp.WorkManager;
import blip.math.Math:sqrt;
import blip.core.sync.Mutex;
import blip.narray.NArray;
import blip.bindings.blas.Types:BlasTypeForType;
import blip.core.Traits: ctfe_i2a;
import dchem.sys.DynVars;
import blip.math.Math: min,max;
import blip.util.RefCount;

/// various levels of duplication
enum PSDupLevel{
    None=0,                 /// does not copy (if possible)
    PSysLevel=1,            /// shallowest level, just duplicates the ParticleSys, but not the internal structures
    DynProperties=2,        /// duplicates dynamic properties
    DynPNullX=4,            /// nullifies x (positions) if DynProperties are duplicated
    DynPNullDx=8,           /// nullifies dx (velocities) if DynProperties are duplicated
    DynPNullMddx=0x10,      /// nullifies mddx (forces) if DynProperties are duplicated
    SysStruct=0x20,         /// duplicates the system structure
    SysStructContents=0x40, /// duplicates the system structure contents
    HiddenVars=0x80,        /// duplicates hidden variables
    DynPNullE=0x100,        /// nullifies the energy if DynProperties are duplicated
    All=PSysLevel|DynProperties|SysStruct|SysStructContents|HiddenVars, /// duplicates all
    EmptyDyn=DynPNullX|DynPNullDx|DynPNullMddx|PSysLevel|DynProperties|HiddenVars, /// duplicates the structure, but uses an empty dynamic info
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
    void addToTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,XType)diffP,DynPVector!(`~t~`,DxType)res);
    /// adds to the vector in res the derivative deriv from the dual tangential space of pos.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order.
    void addFromDualTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,DualDxType)deriv,DynPVector!(`~t~`,XType)res);
    /// transfers from a foreign dual tangential space, this is conceptually equivalent to parallel
    /// transfer, but needs to be valid only if x and pos2.x are close, in particular
    /// a direct transfer res+=derivAtPos2 is always allowed.
    /// uses the dual tangential space as it should be more "uniform"
    void addFromForeignDualTangentialSpace(ParticleSys!(`~t~`)pos1, ParticleSys!(`~t~`)pos2,DynPVector!(`~t~`,DualDxType)derivAtPos2,DynPVector!(`~t~`,DualDxType)res);
    `;
}
/// implements DerivTransferT using the templates that end with T (to work around compiler bugs in implementing interfaces with alias)
char[] derivTransferTMixin(char[]t){
    return `
    void addToTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,XType)diffP,DynPVector!(`~t~`,DxType)res){
        addToTangentialSpaceT!(`~t~`)(pos,diffP,res);
    }
    void addFromDualTangentialSpace(ParticleSys!(`~t~`)pos,DynPVector!(`~t~`,DualDxType)deriv,DynPVector!(`~t~`,XType)res){
        addFromDualTangentialSpaceT!(`~t~`)(pos,deriv,res);
    }
    void addFromForeignDualTangentialSpace(ParticleSys!(`~t~`)pos1,ParticleSys!(`~t~`)pos2,DynPVector!(`~t~`,DualDxType)derivAtPos2,DynPVector!(`~t~`,DualDxType)res){
        addFromForeignDualTangentialSpaceT!(`~t~`)(pos1,pos2,derivAtPos2,res);
    }
    void setOverlap(DynPMatrix!(`~t~`,DxType,DualDxType)m){
        setOverlapT!(`~t~`)(m);
    }
    `;
}

/// interface that implements Real and LowP versions (multiple inheritane of templated interface seems to trigger bugs)
interface DerivTransfer{
    void addToTangentialSpace(ParticleSys!(Real)pos,DynPVector!(Real,XType)diffP,DynPVector!(Real,DxType)res);
    void addToTangentialSpace(ParticleSys!(LowP)pos,DynPVector!(LowP,XType)diffP,DynPVector!(LowP,DxType)res);
    void addFromDualTangentialSpace(ParticleSys!(Real)pos,DynPVector!(Real,DualDxType)deriv,DynPVector!(Real,XType)res);
    void addFromDualTangentialSpace(ParticleSys!(LowP)pos,DynPVector!(LowP,DualDxType)deriv,DynPVector!(LowP,XType)res);
    void addFromForeignDualTangentialSpace(ParticleSys!(Real)pos,ParticleSys!(Real)pos2,DynPVector!(Real,DualDxType)derivAtPos2,DynPVector!(Real,DualDxType)res);
    void addFromForeignDualTangentialSpace(ParticleSys!(LowP)pos,ParticleSys!(LowP)pos2,DynPVector!(LowP,DualDxType)derivAtPos2,DynPVector!(LowP,DualDxType)res);
    /// should return true if the overlap is different from the identity
    bool specialOverlap();
    /// should set the overlap, matrix should use 2d distribution, at the moment it uses just pos,orient,dof sequence
    void setOverlap(DynPMatrix!(BlasTypeForType!(Real),DxType,DualDxType)m);
    static if(! is(BlasTypeForType!(Real)==BlasTypeForType!(LowP))){
        /// should set the overlap, matrix should use 2d distribution, at the moment it uses just pos,orient,dof sequence
        void setOverlap(DynPMatrix!(BlasTypeForType!(LowP),DxType,DualDxType)m);
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
                scope a1=NArray!(ArrBData,2)([cast(index_type)(T.sizeof*blockSizeDarr),ArrBData.sizeof],
                    [arr.length/((blockSizeDarr==0)?1:blockSizeDarr),nElements*T.sizeof/ArrBData.sizeof],0,arr.basicData,0);
                alias typeof(darr.basicData[0]) DarrBData;
                scope a2=NArray!(DarrBData,2)([cast(index_type)(V.sizeof*blockSizeArr),DarrBData.sizeof],
                    [darr.length/((blockSizeArr==0)?1:blockSizeArr),nDelements*V.sizeof/DarrBData.sizeof],0,darr.basicData,0);
                a2+=a1;
            } else {
                alias typeof(darr.basicData[0]) DarrBData;
                assert(el2DelMap.length>1&&(el2DelMap.length-2)%3==0,"unexpected el2DelMap length");
                index_type resEl=blockSizeArr-nElements;
                index_type resDel=blockSizeDarr-nDelements;
                size_t nSeg=el2DelMap.length/3;
                size_t nPart=arr.length/nElements;
                auto p=arr.ptr;
                auto dp=darr.ptr;
                for (size_t iPart=nPart;iPart!=0;--iPart){
                    for(size_t iseg=0;iseg!=nSeg;++iseg){
                        p+=el2DelMap[iseg*3];
                        dp+=el2DelMap[iseg*3+1];
                        for (index_type i=0;i!=el2DelMap[iseg*3+2];--i){
                            for (size_t iEl=0;iEl<V.sizeof/DarrBData.sizeof;++iEl){
                                static if(is(typeof((*dp)+=*p))){
                                    (*dp)+=(*p);
                                } else {
                                    *dp=(*dp)+(*p);
                                }
                                ++dp;
                                ++p;
                            }
                        }
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
                scope a1=NArray!(DarrBData,2)([cast(index_type)(V.sizeof*blockSizeDarr),DarrBData.sizeof],
                    [darr.length/((blockSizeDarr==0)?1:blockSizeDarr),nDelements*(V.sizeof/DarrBData.sizeof)],0,darr.basicData,0);
                alias typeof(arr.basicData[0]) ArrBData;
                scope a2=NArray!(ArrBData,2)([cast(index_type)(T.sizeof*blockSizeArr),ArrBData.sizeof],
                    [arr.length/((blockSizeArr==0)?1:blockSizeArr),nElements*(T.sizeof/ArrBData.sizeof)],0,arr.basicData,0);
                a2+=a1;
            } else {
                alias typeof(arr.basicData[0]) ArrBData;
                assert(el2DelMap.length>1&&(el2DelMap.length-2)%3==0,"unexpected el2DelMap length");
                index_type resEl=blockSizeArr-nElements;
                index_type resDel=blockSizeDarr-nDelements;
                size_t nSeg=el2DelMap.length/3;
                size_t nPart=arr.length/nElements;
                auto p=arr.ptr;
                auto dp=darr.ptr;
                for (size_t iPart=nPart;iPart!=0;--iPart){
                    for(size_t iseg=0;iseg!=nSeg;++iseg){
                        p+=el2DelMap[iseg*3];
                        dp+=el2DelMap[iseg*3+1];
                        for (index_type i=0;i!=el2DelMap[iseg*3+2];--i){
                            for (size_t iEl=0;iEl<T.sizeof/ArrBData.sizeof;++iEl){
                                static if(is(typeof((*dp)+=*p))){
                                    (*p)+=(*dp);
                                } else {
                                    (*p)=(*p)+(*dp);
                                }
                                ++dp;
                                ++p;
                            }
                        }
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
                        // skip what there is to skip
                        index_type elSkipNow=0;
                        if(skipElPos<skipElements.length && elPos>=skipElements[skipElPos]){
                            elSkipNow=skipElements[skipElPos+1]-skipElements[skipElPos];
                            skipElPos+=2;
                        }
                        index_type delSkipNow=0;
                        if(skipDelPos<skipDelements.length && delPos>=skipDelements[skipDelPos]){
                            delSkipNow=skipElements[skipDelPos+1]-skipDelements[skipDelPos];
                            skipDelPos+=2;
                        }
                        elPos+=elSkipNow;
                        delPos+=delSkipNow;
                        
                        index_type nMax=nElements-elPos;
                        if(skipElPos<skipElements.length){
                            if (elPos<skipElements[skipElPos]){
                                nMax=skipElements[skipElPos]-elPos;
                            } else {
                                assert(0,"unexpected elPos>=skipElPos[skipElPos]"); // non merged/invalid skips?
                            }
                        }
                        if (skipDelPos<skipDelements.length){
                            if (delPos<skipDelements[skipDelPos]){
                                nMax=min(nMax,skipDelements[skipDelPos]);
                            } else {
                                assert(0,"unexpected delPos>=skipDelPos[skipDelPos]"); // non merged/invalid skips?
                            }
                        }
                        
                        if (resPos>=el2DelMap.length){
                            el2DelMap~=[elSkipNow,delSkipNow,nMax];
                        } else {
                            assert(el2DelMap.length>resPos+2,"internal error el2DelMap");
                            el2DelMap[resPos]=elSkipNow;
                            el2DelMap[resPos+1]=delSkipNow;
                            el2DelMap[resPos+2]=nMax;
                        }
                        elPos+=nMax;
                        delPos+=nMax;
                        resPos+=3;
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
        foreach (pool;[&p.dynVars.dVarStruct.xGroup]){
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
        foreach (pool;[&p.dynVars.dVarStruct.dxGroup,&p.dynVars.dVarStruct.dualDxGroup]){
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
    void addToTangentialSpaceT(T)(ParticleSys!(T)pos,DynPVector!(T,XType)diffP,DynPVector!(T,DxType)res){
        if (pKind in res.pos.kRange && pKind in diffP.pos.kRange){
            posEls.addArrayToDArray(diffP.pos[pKind],res.pos[pKind],
                diffP.pos.arrayStruct.kindDim(pKind),res.pos.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.orient.kRange && pKind in diffP.orient.kRange){
            orientEls.addArrayToDArray(diffP.orient[pKind],res.orient[pKind],
                diffP.orient.arrayStruct.kindDim(pKind),res.orient.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.dof.kRange && pKind in diffP.dof.kRange){
            dofEls.addArrayToDArray(diffP.dof[pKind],res.dof[pKind],
                diffP.dof.arrayStruct.kindDim(pKind),res.dof.arrayStruct.kindDim(pKind));
        }
        foreach(obj;transferHandlers){
            obj.addToTangentialSpace(pos,diffP,res);
        }
    }
    /// adds to the vector in res the derivative deriv from the tangential space of pos.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order.
    void addFromDualTangentialSpaceT(T)(ParticleSys!(T)pos,DynPVector!(T,DualDxType)deriv,DynPVector!(T,XType)res){
        if (pKind in res.pos.kRange && pKind in deriv.pos.kRange){
            posEls.addDArrayToArray(deriv.pos[pKind],res.pos[pKind],
                deriv.pos.arrayStruct.kindDim(pKind),res.pos.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.orient.kRange && pKind in deriv.orient.kRange){
            orientEls.addDArrayToArray(deriv.orient[pKind],res.orient[pKind],
                deriv.orient.arrayStruct.kindDim(pKind),res.orient.arrayStruct.kindDim(pKind));
        }
        if (pKind in res.dof.kRange && pKind in deriv.dof.kRange){
            dofEls.addDArrayToArray(deriv.dof[pKind],res.dof[pKind],
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
        DynPVector!(T,DualDxType)derivAtPos2,DynPVector!(T,DualDxType)res){
        res.opBypax(derivAtPos2,1,1);
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
    mixin(serializeSome("dchem.SysKind",""));
}

/// represent the structure of a system of particles
class SysStruct: CopiableObjectI,Serializable
{
    char[] name; /// name of the structure
    SubMapping fullSystem; /// LocalIndex is this system struct
    SubMapping externalOrder; /// PIndex is the external one, LocalIndex maps to the PIndex of the current system
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
        assert(particleKinds.contiguous(),"available only for contiguous arrays");
        return particleKinds.support();
    }
    SysStruct dup(PSDupLevel l){
        if ((l&PSDupLevel.SysStructContents)!=0) {
            return new SysStruct(name, fullSystem, externalOrder, levels.dup, particlesStruct, particles.dup,
                superParticle.dup, subParticlesStruct, subParticles.dup, particleKindsStruct, particleKinds.dup);
        } else if ((l&PSDupLevel.SysStruct)!=0){
            return new SysStruct(name, fullSystem, externalOrder, levels, particlesStruct, particles,
                superParticle, subParticlesStruct, subParticles, particleKindsStruct, particleKinds);
        }
        return this;
    }
    SysStruct dup(){
        return dup(PSDupLevel.SysStruct);
    }
    SysStruct deepdup(){
        return dup(PSDupLevel.SysStructContents);
    }
    void reallocDynVarsStructs(T)(ref DynamicsVars!(T) dV){
        dV.giveBack();
        if (dV.dVarStruct !is null) dV.dVarStruct.flush(); // call stopCaching ?
        dV.dVarStruct=new DynamicsVarsStruct!(T)(fullSystem,KindRange(levels[0].kStart,levels[$-1].kEnd));
    }
    mixin(serializeSome("dchem.sys.SysStruct",
        `name: the name of the system
        fullSystem: sub mapping to the whole system
        externalOrder: order of the particles in the external files
        levels: kind ranges of the various levels
        particles: particle indexes
        superParticle: super particle, i.e. molecule for example
        subParticles: the subparticles of each particle
        particleKinds: particle kinds`));
    mixin printOut!();
}

/// interface for hidden variables
interface HiddenVars:Serializable{
    /// duplicates the hidden variables
    HiddenVars dup();
    /// performs a copy of the hidden variables
    HiddenVars emptyCopy();
    /// copies contents of hidden variables (should work even if the argument is null)
    void opSliceAssign(HiddenVars);
    /// updates the hidden variables to a new particle system
    void updateTo(ParticleSys!(Real)newP);
    /// updates the hidden variables to a new particle system
    void updateTo(ParticleSys!(LowP)newP);
}

/// a pool for commonly used particle properties
struct CommonParticleProperties{
    SegArrMemMap!(LowP) poolLowPParticleProperty;
    SegArrMemMap!(Real) poolRealParticleProperty;
    SegArrMemMap!(int)  poolIntParticleProperty;
    SegArrMemMap!(long) poolLongParticleProperty;
    
    void reallocWithParticleStruct(SegmentedArrayStruct pStruct){
        if (poolLowPParticleProperty!is null) poolLowPParticleProperty.rmUser();
        if (poolRealParticleProperty!is null) poolRealParticleProperty.rmUser();
        if (poolIntParticleProperty!is null)  poolIntParticleProperty.rmUser();
        if (poolLongParticleProperty!is null) poolLongParticleProperty.rmUser();
        poolLowPParticleProperty=new SegArrMemMap!(LowP)(pStruct);
        poolRealParticleProperty=new SegArrMemMap!(Real)(pStruct);
        poolIntParticleProperty =new SegArrMemMap!(int )(pStruct);
        poolLongParticleProperty=new SegArrMemMap!(long)(pStruct);
        poolIntParticleProperty.consolidate(poolLowPParticleProperty);
        poolIntParticleProperty.consolidate(poolRealParticleProperty);
        poolLongParticleProperty.consolidate(poolLowPParticleProperty);
        poolLongParticleProperty.consolidate(poolRealParticleProperty);
    }
    SegmentedArray!(T) particlePropertyT(T)(){
        static if (is(T==LowP)){
            return poolLowPParticleProperty.newArray();
        } else static if (is(T==Real)){
            return poolRealParticleProperty.newArray();
        } else static if (is(T==int)){
            return poolIntParticleProperty.newArray();
        } else static if (is(T==long)){
            return poolLongParticleProperty.newArray();
        } else {
            static assert(0,T.stringof~" not supported");
        }
    }
    CommonParticleProperties dup(){
        if (poolLowPParticleProperty!is null) poolLowPParticleProperty.addUser();
        if (poolRealParticleProperty!is null) poolRealParticleProperty.addUser();
        if (poolIntParticleProperty !is null) poolIntParticleProperty .addUser();
        if (poolLongParticleProperty!is null) poolLongParticleProperty.addUser();
        return *this;
    }
    void giveBack(){
        if (poolLowPParticleProperty!is null) poolLowPParticleProperty.rmUser();
        if (poolRealParticleProperty!is null) poolRealParticleProperty.rmUser();
        if (poolIntParticleProperty !is null) poolIntParticleProperty.rmUser();
        if (poolLongParticleProperty!is null) poolLongParticleProperty.rmUser();
    }
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
    
    CommonParticleProperties particlePropertiesPools; /// pool for simple scalar properties
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
                    if (overlap is null || overlap.rowGroup.posStruct !is pSys.dynVars.dVarStruct.dualDxGroup.posStruct
                        || overlap.rowGroup.orientStruct !is pSys.dynVars.dVarStruct.dualDxGroup.orientStruct
                        || overlap.rowGroup.dofStruct !is pSys.dynVars.dVarStruct.dualDxGroup.dofStruct
                        || overlap.colGroup.posStruct !is pSys.dynVars.dVarStruct.dxGroup.posStruct
                        || overlap.colGroup.orientStruct !is pSys.dynVars.dVarStruct.dxGroup.orientStruct
                        || overlap.colGroup.dofStruct !is pSys.dynVars.dVarStruct.dxGroup.dofStruct)
                    {
                        auto rowGroup=pSys.dynVars.dVarStruct.dxGroup;
                        auto colGroup=pSys.dynVars.dVarStruct.dualDxGroup;
                        static if (! is(dtype==dtypeBlas)){
                            auto rGroup=new DynPVectStruct!(dtypeBlas)("blasDxGroup",rowGroup.posStruct,
                                rowGroup.orientStruct,rowGroup.dofStruct);
                            auto cGroup=new DynPVectStruct!(dtypeBlas)("blasDualDxGroup",colGroup.posStruct,
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
                    if (overlapInv is null || overlapInv.rowGroup!is pSys.dynVars.dVarStruct.dxGroup 
                        || overlapInv.colGroup!is pSys.dynVars.dVarStruct.dualDxGroup){
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
        CommonParticleProperties particlePropertiesPools,HiddenVars hVars=null){
        this.iteration=iter;
        this.name=name;
        this.sysStruct=sysStruct;
        this.particlePropertiesPools=particlePropertiesPools;
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
    }
    
    ParticleSys!(V) dupT(V=T)(PSDupLevel l){
        auto dVar=dynVars;
        if ((l&PSDupLevel.DynPNullX))     dVar.x.clear();
        if ((l&PSDupLevel.DynPNullDx))    dVar.dx.clear();
        if ((l&PSDupLevel.DynPNullMddx))  dVar.mddx.clear();
        if ((l&PSDupLevel.DynPNullE))     dVar.potentialEnergy=Real.init;
        static if(is(T==V)){
            if ((l&PSDupLevel.DynProperties)==0){
                return this;
            }
            auto newDynVars=(((l&PSDupLevel.DynProperties)!=0)?dVar.dupT!(V)():dVar);
        } else {
            auto newDynVars=dVar.dupT!(V)();
        }
        return new ParticleSys!(V)(iteration,name,sysStruct.dup(l),null,newDynVars,
            particlePropertiesPools.dup(),
            ((hVars!is null && ((l&PSDupLevel.HiddenVars)!=0))?hVars.dup:hVars));
    }
    
    ParticleSys dup(PSDupLevel l){
        return dupT!(T)(l);
    }
    
    ParticleSys dup(){
        return dupT!(T)(PSDupLevel.DynProperties);
    }
    
    ParticleSys deepdup(){
        return dup(PSDupLevel.All);
    }
    /// reallocates array structs
    void reallocStructs(){
        sysStruct.reallocDynVarsStructs(dynVars);
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
    }
    /// system structure changed (particle added/removed, kinds added/removed)
    /// the segmented array structs should be initialized, and modifiable.
    /// positions,... are not yet available
    void sysStructChanged(){
        this.particlePropertiesPools.reallocWithParticleStruct(sysStruct.particlesStruct);
        derivOverlap.clear();
        this.reallocStructs();
        foreach(pKind;sysStruct.particleKinds.sDataLoop){
            pKind.sysStructChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("sysStructChanged",Variant(this));
        dynVars.dVarStruct.freezeStructs();
    }
    /// position of particles changed, position,... are valid
    void positionsChanged(){
        derivOverlap.invalidate();
        foreach(pKind;sysStruct.particleKinds.sDataLoop){
            pKind.positionsChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("positionsChanged",Variant(this));
    }
    /// cell changed
    void cellChanged(){
        foreach(pKind;sysStruct.particleKinds.sDataLoop){
            pKind.cellChanged(Variant(this));
        }
        if (nCenter!is null)
            nCenter.notify("cellChanged",Variant(this));
    }
    /// copy op
    void opSliceAssign(V)(V p){
        static if(is(V U:ParticleSys!(U))){
            if (p is this) return;
            iteration=p.iteration;
            sysStruct=p.sysStruct;
            dynVars[]=p.dynVars;
            if (hVars!is null) hVars[]=p.hVars;
        } else static if (is(typeof(p.copyTo(this)))){
            p.copyTo(this);
        } else {
            static assert(0,"no assignment possible from type "~V.stringof~" to ParticleSys!("~T.stringof~")");
        }
    }
    
    /// returns a vector in the tangential space (derivative space) equivalent (to first order)
    /// to the one that goes from x-0.5*scale*diffP to x+0.5*scale*diffP
    void addToTSpace(T)(DynPVector!(T,XType)diffP,DynPVector!(T,DxType)res){
        foreach(k;diffP.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addToTangentialSpace(this,diffP,res);
        }
    }
    /// adds to the vector in res the derivative deriv.
    /// This is conceptually equivalent to geodesic following, but it doesn't have to be exactly that,
    /// it just have to be exact at the first order in particular
    /// r[]=0; r2[]=0; addFromDualTangentialSpace(deriv,r); r2[]=r-x; addToTangentialSpace(r2,d);
    /// implies that d and deriv might differ, but should become equal for small derivatives
    /// (in first order).
    void addFromDualTSpace(T)(DynPVector!(T,DualDxType)deriv,DynPVector!(T,XType)res){
        foreach(k;deriv.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addFromDualTangentialSpace(this,deriv,res);
        }
    }
    /// transfers from a foreign tangential space, this is conceptually equivalent to parallel
    /// transfer, but needs to be valid only if x and pos2.x are close, in particular
    /// a direct transfer res+=derivAtPos2 is always allowed.
    void addFromForeignDualTSpace(T)(DynPVector!(T,XType)pos2,DynPVector!(T,DualDxType)derivAtPos2,DynPVector!(T,DualDxType)res){
        foreach(k;derivAtPos2.gRange.intersect(res.gRange).pLoop()){
            sysStruct.particleKinds[k][0].addFromForeignDualTangentialSpace(this,pos2,derivAtPos2,res);
        }
    }
    /// returns either an up to date overlap of the derivatives or null (which means that it is the identity)
    DynPMatrix!(dtypeBlas,DxType,DualDxType) maybeDerivOverlap(){
        if (!derivOverlap.specialOverlap) return null;
        if (derivOverlap.overlapUpToDate) return derivOverlap.overlap;
        derivOverlap.updateOverlap(this);
        if (!derivOverlap.specialOverlap) return null;
        return derivOverlap.overlap;
    }
    /// returns either an up to date inverse of the derivatives overlap or null (which means that it is the identity)
    DynPMatrix!(dtypeBlas,DualDxType,DxType) maybeDerivOverlapInv(){
        if (!derivOverlap.specialOverlap) return null;
        if (derivOverlap.overlapInvUpToDate) return derivOverlap.overlapInv;
        derivOverlap.updateOverlapInv(this);
        if (!derivOverlap.specialOverlap) return null;
        return derivOverlap.overlapInv;
    }
    /// returns true if the overlap of derivatives is not the unit matrix
    bool specialDerivOverlap(){
        bool specOv=false;
        foreach(pKind;sysStruct.particleKinds.sDataLoop){
            if (pKind.specialOverlap()){
                return true;
            }
        }
        return false;
    }
    /// sets the overlap matrix (should be initialized with 0 at the beginnng)
    void setDerivOverlap(DynPMatrix!(dtypeBlas,1,2) overlap){
        foreach(pKind;sysStruct.particleKinds.sDataLoop){
            pKind.setOverlap(overlap);
        }
    }
    
    /// goes to the dual tangential space (multiplication with S^-1) and adds the result to dualDeriv
    void toDualTSpace(DynPVector!(T,DxType)deriv,DynPVector!(T,DualDxType)dualDeriv,T scaleRes=1,T scaleDualDeriv=0){
        auto overlapInv=maybeDerivOverlapInv();
        if (overlapInv!is null){
            overlapInv.matVectMult(deriv,dualDeriv,scaleRes,scaleDualDeriv);
        } else { // assumes two groups ar actually equal
            dualDeriv.opBypax(deriv.toGroup!(2)(),scaleRes,scaleDualDeriv);
        }
    }
    /// comes back from the dual tangential space (multiplication with S)
    void fromDualTSpace(DynPVector!(T,DualDxType)dualDeriv,DynPVector!(T,DxType)deriv,T scaleRes=1,T scaleDeriv=0){
        auto overlap=maybeDerivOverlap();
        if (overlap!is null){
            overlap.matVectMult(dualDeriv,deriv,scaleRes,scaleDeriv);
        } else {
            deriv.opBypax(dualDeriv.toGroup!(1)(),scaleRes,scaleDeriv);
        }
    }
    /// projects into the correct movement space (back & forth from the direct space).
    /// useful to get rid of the components in the null space
    T projectInDualTSpace(DynPVector!(T,DualDxType)dualDeriv){
        T res;
        auto overlap=maybeDerivOverlap();
        if (overlap!is null){
            scope newDirDirect=dynVars.dVarStruct.emptyDx();
            fromDualTSpace(dualDeriv,newDirDirect);
            toDualTSpace(newDirDirect,dualDeriv);
            res=dotInTSpace(newDirDirect,dualDeriv);
            newDirDirect.giveBack();
        } else {
            res=dualDeriv.norm2();
        }
        return res;
    }
    /// projects into the correct movement space (back & forth from the dual space).
    /// useful to get rid of the components in the null space
    /// returns the norm of the vector
    T projectInTSpace(DynPVector!(T,DxType)deriv){
        T res;
        auto overlap=maybeDerivOverlap();
        if (overlap!is null){
            auto dualDeriv=dynVars.dVarStruct.emptyDualDx();
            toDualTSpace(deriv,dualDeriv);
            fromDualTSpace(dualDeriv,deriv);
            res=dotInTSpace(deriv,dualDeriv);
            dualDeriv.giveBack();
        } else {
            res=deriv.norm2();
        }
        return res;
    }
    /// perform scalar product in the tangential (derivative) space
    T dotInTSpace(DynPVector!(T,DxType)v1,DynPVector!(T,DualDxType)v2){
        assert(v1.pos is null || v2.pos is null || v1.pos.arrayStruct==v2.pos.arrayStruct,"redistribution of vectors not implemented");
        assert(v1.orient is null || v2.orient is null || v1.orient.arrayStruct==v2.orient.arrayStruct,"redistribution of vectors not implemented");
        assert(v1.dof is null || v2.dof is null || v1.dof.arrayStruct==v2.dof.arrayStruct,"redistribution of vectors not implemented");
        return v1.opDot(v2.toGroup!(1)());
    }
    /// updates the hidden vars to the current position
    void updateHVars(){
        if (hVars!is null){
            static if (is(T==Real)||is(T==LowP)){
                hVars.updateTo(this);
            } else {
                // costly!!!
                scope tmpPSys=dupT!(Real)(PSDupLevel.DynProperties);
                tmpPSys.updateHVars();
            }
        }
    }
    
    void release0(){
        dynVars.giveBack();
        particlePropertiesPools.giveBack();
    }
    mixin RefCountMixin!();
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

