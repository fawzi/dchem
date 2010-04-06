module dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import blip.container.BulkArray;
import dchem.sys.SubMapping;
import dchem.Common;
import blip.narray.NArray;
import blip.sync.Atomic;
import blip.t.core.Traits;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import blip.t.core.Variant;
import blip.parallel.smp.WorkManager;
import blip.container.AtomicSLink;
import blip.io.Console;

enum ParaFlags{
    FullPara,
    KindPara,
    Sequential,
}

/// structure of a segmented array, kindStarts is valid only after freezing
final class SegmentedArrayStruct{
    enum Flags{
        None=0,   /// no special flags
        Frozen=1, /// kind dimensions cannot be changed anymore (this is set when an array is created)
        Direct=2, /// =submapping.mappingKind & MappingKind.Direct, worth caching?
        Min1=4,   /// at least one value per kind should be stored
    }
    char[] name; /// a name (for debugging)
    SubMapping submapping;
    KindRange  kRange;
    index_type[]   _kindDims;
    index_type[]   kindStarts;
    Flags flags;
    mixin(serializeSome("dchem.sys.SegmentedArrayStruct","submapping|kRange|_kindDims|kindStarts"));
    mixin printOut!();
    /// allocates a new SegmentedArray with the given kind dimensions
    /// if the Min1 flag is set (default) then at least one value per kind is stored 
    this(char[]name,SubMapping submapping,KindRange kRange,index_type[]kindDims=null,Flags f=Flags.Min1){
        assert(submapping !is null,"submapping needed to allocate");
        this.name=name;
        this.submapping = submapping;
        this.kRange     = kRange    ;
        this.flags      = f;
        this.flags=(f& ~Flags.Direct)|(((submapping.mappingKind & MappingKind.Direct)!=0)?
                Flags.Direct:Flags.None);
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        if (kindDims.length==0 && nkinds!=0){
            kindDims=new index_type[](nkinds);
            kindDims[]=0; // should be the default
        }
        this._kindDims   = kindDims  ;
        assert(_kindDims.length==nkinds);
        assert(kRange in submapping.lKRange,"submapping is smaller than current kRange");
        if ((flags & Flags.Frozen)!=0){
            recalculateStarts();
        } else {
            kindStarts[]=index_type.max;
        }
    }
    /// ditto
    this(char[] name,SubMapping submapping, KindRange kRange, NArray!(index_type,1) kindDims,Flags f=Flags.Min1){
        if (kindDims is null) {
            this(name,submapping,kRange,cast(index_type[])null,f);
        } else {
            index_type[] kDims=new index_type[kindDims.shape[0]];
            for (size_t i=0;i<kindDims.shape[0];++i){
                kDims[i]=kindDims[i];
            }
            this(name,submapping,kRange,kDims,f);
        }
    }
    // internal for serialization
    this(){ }
    /// makes the structure non modifiable (done before creating data arrays)
    typeof(this) freeze(){
        if ((flags & Flags.Frozen)==0){
            synchronized(this){
                if ((flags & Flags.Frozen)==0){
                    recalculateStarts();
                    flags|=Flags.Frozen;
                }
            }
        }
        return this;
    }
    /// recalculates the starts of the array
    void recalculateStarts(){
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        kindStarts.length=nkinds+1;
        
        kindStarts[0]=0;
        bool min1=((flags & Flags.Min1)!=0);
        auto kindShift=cast(size_t)(kRange.kStart-submapping.lKRange.kStart);
        for (size_t i=0;i<nkinds;++i){
            if (min1 && _kindDims[i]==0){
                kindStarts[i+1]=kindStarts[i]+1; // always keep a "kind owned" value
            } else {
                kindStarts[i+1]=kindStarts[i]
                    +_kindDims[i]*(submapping.kindStarts[i+1+kindShift]-submapping.kindStarts[i+kindShift]);
            }
        }
    }
    /// equality, add also equivalence and inclusion??
    equals_t opEquals(Object o2){
        SegmentedArrayStruct s2=cast(SegmentedArrayStruct)o2;
        if (s2!is null){
            return this is s2 || (kRange==s2.kRange && kindStarts==s2.kindStarts
                && submapping==s2.submapping && flags==s2.flags);
        }
        return false;
    }
    /// return the dimension of kind i
    index_type kindDim(T)(T k){
        assert(cast(KindIdx)k in kRange,"kind out of range");
        return _kindDims[cast(size_t)(k-kRange.kStart)];
    }
    /// sets the dimension of kind i, dangerous when called concurrently
    void setKindDim(KindIdx k,index_type val){
        if((flags & Flags.Frozen)==0)
            throw new Exception("index assign on frozen struct",__FILE__,__LINE__);
        assert(k in kRange,"kind out of range");
        auto i=cast(size_t)(k-kRange.kStart);
        return _kindDims[i]=val;
    }
    /// adds the the dimension, returns the original value
    index_type addToKindDim(KindIdx k,index_type val){
        assert(k in kRange,"kind out of range");
        auto i=cast(size_t)(k-kRange.kStart);
        if((flags & Flags.Frozen)==0)
            throw new Exception("index assign on frozen struct",__FILE__,__LINE__);
        return atomicAdd(_kindDims[cast(size_t)i],val);
    }
    /// returns a non frozen copy
    typeof(this) dup(){
        return new SegmentedArrayStruct(name~"Dup",submapping,kRange,_kindDims.dup,flags & (~Flags.Frozen));
    }
}


/// segmented (level/kinds) array (with kind uniform dimension)
final class SegmentedArray(T){
    SegmentedArrayStruct arrayStruct;
    KindRange  kRange;
    index_type[]   kindStarts;
    BulkArray!(T) _data;
    bool direct;
    alias T dtype;
    static size_t defaultOptimalBlockSize=32*1024/T.sizeof;
    
    mixin(serializeSome("dchem.sys.SegmentedArray","kRange|kindStarts|_data"));
    mixin printOut!();
    
    BulkArray!(T) data(){
        return _data[kindStarts[0],kindStarts[$-1]];
    }

    // internal for serialization
    this(){ }
    /// allocates a new SegmentedArray with the given kind dimensions
    /// min1 
    this(SegmentedArrayStruct arrayStruct, BulkArray!(T)data=BulkArray!(T).dummy,
        KindRange kRange=KindRange.all,index_type[] kindStarts=null)
    in{
        auto myKRange=kRange;
        if (kRange.kEnd==KindIdx.init) myKRange=arrayStruct.kRange;
        if (kindStarts!is null){
            auto nkinds=cast(size_t)(myKRange.kEnd-myKRange.kStart);
            assert(kindStarts.length==nkinds+1,"kindStarts has wrong size");
            auto kindShift=cast(size_t)(myKRange.kStart-arrayStruct.kRange.kStart);
            
            for (size_t i=0;i<=nkinds;++i){
                assert(kindStarts[i]==kindStarts[0]+arrayStruct.kindStarts[kindShift+i]-arrayStruct.kindStarts[kindShift],
                    "invalid kindStarts");
            }
        }
    } body {
        assert(arrayStruct !is null,"arrayStruct must be valid");
        this.arrayStruct=arrayStruct.freeze;
        this.kRange=kRange;
        if (kRange.kEnd==KindIdx.init){
            assert(kRange.kStart==0);
            this.kRange=arrayStruct.kRange;
        }
        assert(this.kRange in arrayStruct.kRange);
        this.kindStarts=kindStarts;
        if (this.kindStarts is null){
            auto nkinds=cast(size_t)(this.kRange.kEnd-this.kRange.kStart);
            this.kindStarts=new index_type[cast(size_t)(nkinds+1)];
            auto kindShift=cast(size_t)(this.kRange.kStart-arrayStruct.kRange.kStart);
            for (size_t i=0;i<=nkinds;++i){
                this.kindStarts[i]=arrayStruct.kindStarts[kindShift+i]-arrayStruct.kindStarts[kindShift];
            }
        }
        _data=data;
        if (BulkArrayIsDummy(data)){
            _data=BulkArray!(T)(this.kindStarts[this.kindStarts.length-1]-this.kindStarts[0]);
        }
        direct=(arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0;
    }
    /// array of elements in the given range
    SegmentedArray!(T) opIndex(KindRange kr){
        KindRange krCommon=kRange.intersect(kr);
        size_t startLK=cast(size_t)(krCommon.kStart-kRange.kStart);
        size_t endLK=cast(size_t)(krCommon.kEnd-kRange.kStart+1);
        return new SegmentedArray(arrayStruct,data,krCommon,kindStarts[startLK..endLK]);
    }
    /// array of elements for *local* kind k
    BulkArray!(T) opIndex(KindIdx k){
        auto res=data[kindStarts[cast(size_t)k],kindStarts[cast(size_t)k+1]];
        // res.blockSize=arrayStruct.kindDim[cast(size_t)k];
        return res;
    }
    /// gets the particle using the local particle numbering, but that might be outside the kind range
    T[] getMaybeInRange(LocalPIndex p){
        auto kindIdx=p.kind();
        if (kindIdx in kRange){
            auto startIdx=kindStarts[cast(size_t)kindIdx]+p.particle();
            return _data.getSlice(startIdx,startIdx+arrayStruct.kindDim(kindIdx));
        }
        return null;
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    T[] opIndex(LocalPIndex p){
        auto kindIdx=p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
        auto dim=arrayStruct.kindDim(kindIdx);
        auto startIdx=kindStarts[cast(size_t)(kindIdx-kRange.kStart)]+p.particle()*dim;
        return _data.getSlice(startIdx,startIdx+((dim!=0)?dim:
            (arrayStruct.flags&SegmentedArrayStruct.Flags.Min1!=0)?1:0));
    }
    /// array of elements for a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T[] opIndex(PIndex p){
        if ((arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0){
            auto kindIdx=p.kind();
            auto pos=p.particle();
            assert(kindIdx in kRange,"kind out of range");
            auto dim=arrayStruct.kindDim(kindIdx);
            auto startIdx=kindStarts[cast(size_t)kindIdx-kRange.kStart]+p.particle()*dim;
            assert(startIdx<kindStarts[kindIdx-kRange.kStart+1] || dim==0,"particle index out of bounds");
            return _data.getSlice(startIdx,startIdx+((dim!=0)?dim:
                ((arrayStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0)?1:0));
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            auto dim=arrayStruct.kindDim(kindIdx);
            auto startIdx=kindStarts[kindIdx-kRange.kStart]+l.particle()*dim;
            assert(startIdx<kindStarts[kindIdx-kRange.kStart+1] || dim==0,"particle index out of bounds");
            return _data.getSlice(startIdx,startIdx+((dim!=0)?dim:
                ((arrayStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0)?1:0));
        }
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T *ptrI(LocalPIndex p,index_type i){
        auto kindIdx=p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
        assert(i>=0 && (i<arrayStruct.kindDim(kindIdx)||
                (i==0&&(arrayStruct.flags & SegmentedArrayStruct.Flags.Min1)!=0)),"index i out of bounds");
        auto kPos=cast(size_t)(kindIdx-kRange.kStart);
        auto dim=arrayStruct.kindDim(kindIdx);
        auto startIdx=kindStarts[kPos]+p.particle()*dim;
        assert(startIdx+dim<=kindStarts[kPos+1],"p.particle out of bounds");
        return _data.ptrI(startIdx+i);
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T opIndex(LocalPIndex p,index_type i){
        return *ptrI(p,i);
    }
    /// sets the value for element i of the particle using the local particle numbering
    /// (has to be in the kind range)
    void opIndexAssign(T val,LocalPIndex p,index_type i){
        static if (is(typeof(val.opSliceAssign(val)))){
            (*ptrI(p,i))[]=val;
        } else {
            *ptrI(p,i)=val;
        }
    }
    /// address of element i of a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T *ptrI(PIndex p,index_type i){
        if (direct){
            auto kindIdx=p.kind();
            assert(kindIdx in kRange,"kind out of range");
            auto dim=arrayStruct.kindDim(kindIdx);
            assert(i>=0 && (i<dim||(i==0&&(arrayStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0)),"index i out of bounds");
            auto pos=p.particle();
            auto startIdx=kindStarts[kindIdx-kRange.kStart]+p.particle()*dim;
            assert(startIdx+dim<=kindStarts[kindIdx-kRange.kStart+1],"p.particle out of range");
            return _data.ptrI(startIdx+cast(size_t)i);
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            auto dim=arrayStruct.kindDim(kindIdx);
            assert(i>=0 && (i<dim||(i==0&&(arrayStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0)),"index i out of bounds");
            auto startIdx=kindStarts[cast(size_t)(kindIdx-kRange.kStart)]+l.particle()*dim;
            assert(startIdx+dim<=kindStarts[kindIdx-kRange.kStart+1],"p.particle out of range");
            return _data.ptrI(startIdx+i);
        }
    }
    /// element i of a (global) particle p array
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T opIndex(PIndex p,index_type i){
        return *ptrI(p,i);
    }
    /// copies from an array to this
    void opSliceAssign(V)(SegmentedArray!(V) val){
        assert(arrayStruct==val.arrayStruct,"different structs");
        assert(kRange==val.kRange,"different kRanges");
        assert(kindStarts==val.kindStarts,"different kindStarts");
        _data.copyFrom(val._data);
    }
    /// copies from an array to this
    void opSliceAssign()(T val){
        _data[]=val;
    }
    /// copies this array to the given SegmentedArray, tryig to reuse its memory allocations
    void dupTo(V)(SegmentedArray!(V) val){
        val.arrayStruct=arrayStruct;
        val.kRange=kRange;
        val.kindStarts=kindStarts;
        val.direct=direct;
        static if(is(T==V)){
            if (_data.length!=val._data.length){
                val._data.dataOfGuard(_data.dup().guard);
            } else {
                val._data.copyFrom!(V)(_data);
            }
        } else {
            if (_data.length!=val._data.length){
                val._data.dataOfGuard(new blip.container.BulkArray.Guard(_data.length*V.sizeof));
            } else {
                val._data.copyFrom!(T)(_data);
            }
        }
    }
    /// returns a copy of the segmented array
    SegmentedArray!(V) dupT(V=T)(){
        return new SegmentedArray!(V)(arrayStruct,_data.dupT!(V),kRange,kindStarts.dup);
    }
    SegmentedArray dup(){
        return new SegmentedArray(arrayStruct,_data.dup,kRange,kindStarts.dup);
    }
    /// returns a copy of the segmented array
    SegmentedArray deepdup(){
        return new SegmentedArray(arrayStruct,_data.deepdup,kRange,kindStarts.dup);
    }
    
    /// loops on the particles and each element of a multivaluated segmented array.
    /// for consistency always skips kinds wiouth particles (even if Min1)
    struct PLoop(int pFlags){
        size_t optimalBlockSize;
        SegmentedArray array;
        int opApply(int delegate(ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","if (auto res=dlg(*arrayPtr)){ *finalRes=res; return; }","",
            
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayNel;++ii){
                if (auto res=dlg(arrayPtr[ii])){ *finalRes=res; return; }
            }`,""]));
            return result;
        }
        int opApply(int delegate(ref LocalPIndex lIdx,ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref LocalPIndex,ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","if (auto res=dlg(localPIndex,*arrayPtr)){ *finalRes=res; return; }","",
            
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayNel;++ii){
                if (auto res=dlg(localPIndex,arrayPtr[ii])){ *finalRes=res; return; }
            }`,""]));
            return result;
        }
        int opApply(int delegate(ref PIndex pIdx,ref LocalPIndex lIdx,ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref PIndex, ref LocalPIndex, ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","if (auto returnV=dlg(*pIndexPtr,localPIndex,*arrayPtr)){ *finalRes=returnV; return; }","",
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayNel;++ii){
                if (auto returnV=dlg(*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
            }`,""]));
           return result;
        }
        int opApply(int delegate(ref size_t i,ref PIndex pIdx,ref LocalPIndex lIdx,ref T)loopBody){
            int result=0;
            mixin(segArrayMonoLoop(pFlags,"iterContext",["array"],
            "int delegate(ref size_t i,ref PIndex, ref LocalPIndex, ref T) dlg; int* finalRes;","",
            "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtrStart !is null);","",
            ["if ((*finalRes)!=0) return;","size_t ii=0; if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,*arrayPtr)){ *finalRes=returnV; return; }","",
            
            "if ((*finalRes)!=0) return;",`
            for (size_t ii=0;ii<arrayNel;++ii){
                if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
            }`,""]));
           return result;
        }
    }
    /// full parallel loop
    PLoop!(ParaFlags.FullPara) pLoop(size_t optSize=defaultOptimalBlockSize){
        PLoop!(ParaFlags.FullPara) res;
        res.optimalBlockSize=optSize;
        res.array=this;
        return res;
    }
    /// loop parallel only between kinds
    PLoop!(ParaFlags.KindPara) pLoopKinds(size_t optSize=defaultOptimalBlockSize){
        PLoop!(ParaFlags.KindPara) res;
        res.optimalBlockSize=optSize;
        res.array=this;
        return res;
    }
    /// sequential loop
    PLoop!(ParaFlags.Sequential) sLoop(){
        PLoop!(ParaFlags.Sequential) res;
        res.optimalBlockSize=defaultOptimalBlockSize;
        res.array=this;
        return res;
    }
    
    void axpby(V)(SegmentedArray!(V) x,V a=cscalar!(V,1),T b=cscalar!(T,1)){
        auto optimalBlockSize=defaultOptimalBlockSize;
        auto y=this;
        if(b==1){
            if(a==1){
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["y","x"],"","",
                "",`
                visitKind=visitKind&&(!outOfRange);
                assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in axpby");`,"",
                ["","*yPtr += scalar!(T)(*xPtr);","",
                
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    yPtr[ii] += scalar!(T)(xPtr[ii]);
                }`,""]));
            } else {
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["y","x"],
                "V a;","",
                "mainContext.a=a;",`
                visitKind=visitKind&&(!outOfRange);
                assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in axpby");`,"",
                ["","*yPtr += scalar!(T)((*xPtr)*a);","",
                
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    yPtr[ii] += scalar!(T)(xPtr[ii]*a);
                }`,""]));
            }
        } else if(b==0){
            mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["y","x"],
            "V a;","",
            "mainContext.a=a;",`
            visitKind=visitKind&&(!outOfRange);
            assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in axpby");`,"",
            ["","*yPtr = scalar!(T)((*xPtr)*a);","",
            
            "",`
            for (size_t ii=0;ii<xNel;++ii){
                yPtr[ii] = scalar!(T)(xPtr[ii]*a);
            }`,""]));
        } else {
            mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["y","x"],
            "V a;T b;","",
            "mainContext.a=a;mainContext.b=b;",`
            visitKind=visitKind&&(!outOfRange);
            assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in axpby");`,"",
            ["","*yPtr = scalar!(T)((*yPtr)*b+(*xPtr)*a);","",
            
            "",`
            for (size_t ii=0;ii<xNel;++ii){
                yPtr[ii] = scalar!(T)((*yPtr)*b+(*xPtr)*a);
            }`,""]));
        }
    }
    
    static if (is(typeof(T.init*T.init))){
        void opMulAssign()(T scale){
            auto optimalBlockSize=defaultOptimalBlockSize;
            auto x=this;
            static if (is(typeof(function(T t){ t *= t; }))){
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["x"],
                "T scale;","",
                "mainContext.scale=scale;","visitKind=visitKind&&(!outOfRange);","",
                ["","*xPtr *= scale;","",
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    xPtr[ii] *= scale;
                }`,""]));
            } else {
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["x"],
                "T scale;","",
                "mainContext.scale=scale;","visitKind=visitKind&&(!outOfRange);","",
                ["","*xPtr = (*xPtr)*scale;","",
                
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    xPtr[ii] = xPtr[ii]*scale;
                }`,""]));
            }
        }
        void opMulAssign(V)(SegmentedArray!(V) y){
            auto optimalBlockSize=defaultOptimalBlockSize;
            auto x=this;
            static if (is(typeof(function(T t){ t *= t; }))){
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["x","y"],
                "T scale;","",
                "mainContext.scale=scale;",`
                visitKind=visitKind&&(!outOfRange);
                assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in opMulAssign");`,"",
                ["","*xPtr *= *yPtr;","",
                
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    xPtr[ii] *= yPtr[ii];
                }`,""]));
            } else {
                mixin(segArrayMonoLoop(ParaFlags.FullPara,"iterContext",["x","y"],
                "T scale;","",
                "mainContext.scale=scale;",`
                visitKind=visitKind&&(!outOfRange);
                assert((kNew.xNel<=1&&kNew.yNel<=1)||kNew.xNel==kNew.yNel,"variable combination of Nel not implemented in opMulAssign");`,"",
                ["","*xPtr = (*xPtr)*(*yPtr);","",
                
                "",`
                for (size_t ii=0;ii<xNel;++ii){
                    xPtr[ii] = xPtr[ii]*yPtr[ii];
                }`,""]));
            }
        }
    }
    
}

/// inner loop variant (exec* method)
char[] segArrayContextExecStr(int pFlags,char[] iterName, char[][]namesLocal,char[] contextExtra,char[] uniq,char[] execN,
    char[]pVisitLocalStart,char[]pVisit,char[]pVisitLocalEnd){
    char[]visitorStr=`
        void exec`~execN~`(){
            if (context.exception !is null) return;`;
    visitorStr~="\n";
    visitorStr~=pVisitLocalStart;
    visitorStr~="\n";
    if (pFlags==ParaFlags.FullPara){
        visitorStr~=`
            if (end-start>optimalBlockSize*3/2){
                auto mid=cast(ParticleIdx)((end-start)/2);
                if (mid>optimalBlockSize){ // tries to make optimalBlockSize a possible fast path
                    mid=cast(ParticleIdx)(((mid+optimalBlockSize-1)/optimalBlockSize)*optimalBlockSize);
                }
                auto firstHalf=alloc();
                firstHalf.end=start+mid;
                Task("SegArrLoop1`~iterName~`",&firstHalf.exec`~execN~`).appendOnFinish(&firstHalf.giveBack).autorelease.submit; // would submitYield be better (less suspended tasks, but more suspensions)??? 
                auto secondHalf=alloc();
                secondHalf.start=start+mid;`;
        foreach (name;namesLocal){
            visitorStr~=`
                secondHalf.`~name~`PtrStart+=`~name~`Nel*mid;`;
        }
        visitorStr~=`
                Task("SegArrLoop2`~iterName~`",&secondHalf.exec`~execN~`).appendOnFinish(&secondHalf.giveBack).autorelease.submit;
            } else {`;
    } else {
        visitorStr~=`
            {`;
    }
    visitorStr~=`
                try{
                    LocalPIndex localPIndex`~uniq~`=LocalPIndex(kind,start);
                    PIndex * pIndexPtr`~uniq~`=pIndexPtrStart;`;
    foreach (name;namesLocal){
        visitorStr~=`
                    `~name~`.dtype* `~name~`Ptr`~uniq~`=`~name~`PtrStart;`;
    }
    visitorStr~=`
                    for (ParticleIdx index`~uniq~`=start;index`~uniq~`<end;++index`~uniq~`){
                        lIndex=index`~uniq~`;
                        `~pVisit~`
                        ++localPIndex`~uniq~`;
                        ++pIndexPtr`~uniq~`;`;
    foreach (name;namesLocal){
        visitorStr~=`
                        `~name~`Ptr`~uniq~`+=`~name~`Nel;`;
    }
    visitorStr~=`
                    }
                } catch (Exception e){
                    context.exception=e;
                }
            }
            `~pVisitLocalEnd~`
        }`;
    return visitorStr;
}

/// the whole context struct, can have several local visitors (pVisitors), those are
/// always sets of 3 strings: pVisitLocalStart,pVisit,pVisitLocalEnd, and they create
/// each one exec context: exec0,exec1,...
char[] segArrayContextStr(int pFlags,char[] iterName, char[][]namesLocal,char[] contextExtra,
    char[] uniq,char[][] pVisitors=[])
{
    assert(pVisitors.length%3==0,"otherVisitors is supposed to be a multiple of 3");
    char[] visitorStr=`
    struct `~iterName~`{`;
    foreach (name;namesLocal){
        visitorStr~=`
        `~name~`.dtype* `~name~`PtrStart;
        size_t `~name~`Nel;`;
    }
    visitorStr~=`
        size_t optimalBlockSize;
        KindIdx kind;
        ParticleIdx start;
        ParticleIdx lIndex;
        ParticleIdx end;
        PIndex *pIndexPtrStart;
        size_t maxInnerDim;
        Exception exception;
        `~iterName~`* context;
        `~iterName~`* next;`;
    visitorStr~="\n";
    visitorStr~=contextExtra;
    visitorStr~=`
        typeof(this) alloc(){
            auto res=popFrom(context.next);
            if (res is null){
                res=new `~iterName~`;
            }
            *res=*this;
            res.next=this;
            return res;
        }`;
    for(int i=0;i<pVisitors.length;i+=3){
        visitorStr~=segArrayContextExecStr(pFlags,iterName,namesLocal,contextExtra,uniq,ctfe_i2a(i/3),
            pVisitors[i],pVisitors[i+1],pVisitors[i+2]);
    }
    visitorStr~=`
        void giveBack(){
            insertAt(context.next,this);
        }
    }
    `;
    return visitorStr;
}

/// loop on the kinds (main external loop)
char[] segArrayKLoopStr(int pFlags,char[] iterName, char[][]namesLocal,
    char[]startLoop,char[] loopBody,char[]endLoop,char[] uniq="")
{
    char[] visitorStr=`
    {
        `~iterName~` mainContext`~uniq~`;
        auto kEnd`~uniq~`=`~namesLocal[0]~`.kRange.kEnd;
        `~iterName~`* newK`~uniq~`;
        PIndex dummyP`~uniq~`;
        mainContext`~uniq~`.optimalBlockSize=optimalBlockSize`~uniq~`;
        mainContext`~uniq~`.context=&mainContext`~uniq~`;`;
    visitorStr~=startLoop;
    if (pFlags!=ParaFlags.Sequential){
        visitorStr~=`
        Task("segArrayKLoop`~iterName~`",delegate void(){`;
    }
    visitorStr~=`
            for (KindIdx kIdx`~uniq~`=`~namesLocal[0]~`.kRange.kStart;kIdx`~uniq~`<kEnd`~uniq~`;++kIdx`~uniq~`){
                if (mainContext.exception !is null) break;
                bool visitKind`~uniq~`=true;
                bool outOfRange`~uniq~`=false;
                if (newK`~uniq~` is null){
                    newK`~uniq~`=mainContext`~uniq~`.alloc();
                }
                newK`~uniq~`.kind=kIdx`~uniq~`;
                newK`~uniq~`.start=0;
                newK`~uniq~`.end=
                    `~namesLocal[0]~`.arrayStruct.submapping.nLocalParticles(kIdx`~uniq~`);
                auto submap`~uniq~`=`~namesLocal[0]~`.arrayStruct.submapping;
                if (submap`~uniq~`.kindStarts[kIdx`~uniq~`-submap`~uniq~`.lKRange.kStart]<
                    submap`~uniq~`.kindStarts[kIdx`~uniq~`-submap`~uniq~`.lKRange.kStart+1])
                {
                    newK`~uniq~`.pIndexPtrStart=submap`~uniq~`.ptrI(
                        LocalPIndex(newK`~uniq~`.kind,newK`~uniq~`.start));
                } else {
                    newK`~uniq~`.pIndexPtrStart=null;
                }`;
    foreach (i,name;namesLocal){
        visitorStr~=`
                if (kIdx`~uniq~` in `~name~`.kRange){
                    auto ik`~uniq~`=cast(size_t)(kIdx`~uniq~`-`~name~`.kRange.kStart);
                    if (`~name~`.kindStarts[ik`~uniq~`]<`~name~`.kindStarts[ik`~uniq~`+1]){
                        newK`~uniq~`.`~name~`PtrStart=`~name~`._data.ptr+`~name~`.kindStarts[ik`~uniq~`];
                        newK`~uniq~`.`~name~`Nel=`~name~`.arrayStruct.kindDim(kIdx`~uniq~`);
                    } else {
                        outOfRange`~uniq~`=true;
                        newK`~uniq~`.`~name~`PtrStart=null;
                        newK`~uniq~`.`~name~`Nel=0;
                    }
                } else {
                    outOfRange`~uniq~`=true;
                    newK`~uniq~`.`~name~`PtrStart=null;
                    newK`~uniq~`.`~name~`Nel=0;
                }`;
    }
    visitorStr~="\n";
    visitorStr~=loopBody;
    visitorStr~="\n";
    visitorStr~=`
        }`;
    if (pFlags!=ParaFlags.Sequential){
        visitorStr~=`
        }).autorelease.executeNow();`;
    }
    visitorStr~=`
        if (mainContext.exception !is null) throw new Exception("exception in SegmentedArray loop",__FILE__,__LINE__,mainContext.exception);
        auto freeL=mainContext`~uniq~`.next;
        while (freeL!is null){
            auto cNext=freeL.next;
            delete freeL;
            freeL=cNext;
        }
    `;
    visitorStr~=endLoop;
    visitorStr~=`
    }`;
    return visitorStr;
}

/// needs "optimalBlockSize"~uniq
/// public vars in kindVisit: struct iterName, "context"~uniq, "kIdx"~uniq, "outOfRange"~uniq, "visitKind"~uniq
/// public vars in pVisit: name~"Ptr"~uniq, "localPIndex"~uniq, "pIndexPtr"~uniq, "index"~uniq
/// there can be several versions of pVisit, if only one is given then that should be the generic version,
/// if two versions are given then the first is when there is one element per particle or per kind,
/// and the other is the generic, if 3 are given there is also one when they have the same size or 1
/// (i.e a single loop with increments 0 and 1 would work)
char[] segArrayMonoLoop(int pFlags,char[] iterName, char[][]namesLocal,
    char[] contextExtra,char[] uniq,char[]startLoop,char[]kindVisit,char[]endLoop,
    char[][]pVisitors)
{
    char[] res="{\n";
    res~=segArrayContextStr(pFlags,iterName, namesLocal, contextExtra,uniq,
        pVisitors);
    char[] kVisit=kindVisit~`
            if (visitKind`~uniq~`){
                newK.maxInnerDim=1;
                bool multiDim`~uniq~`=false;`;
    foreach (n;namesLocal){
                kVisit~=`
                if (newK`~uniq~`.`~n~`Nel!=newK.maxInnerDim && newK`~uniq~`.`~n~`Nel!=0 && newK`~uniq~`.`~n~`Nel!=1){
                    if (newK.maxInnerDim!=1){
                        multiDim`~uniq~`=true;
                    }
                    if (newK.maxInnerDim<newK`~uniq~`.`~n~`Nel) newK.maxInnerDim=newK`~uniq~`.`~n~`Nel;
                }`;
    }
    if (pVisitors.length==3) {
        if (pFlags==ParaFlags.Sequential){
            kVisit~=`
                newK`~uniq~`.exec0();`;
        } else {
            kVisit~=`
                Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                    .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                newK`~uniq~`=null;`;
        }
    } else if (pVisitors.length==6){
        if (pFlags==ParaFlags.Sequential){
            kVisit~=`
                if (newK.maxInnerDim==1){
                    newK`~uniq~`.exec0();
                } else {
                    newK`~uniq~`.exec1();
                }
                newK`~uniq~`=null;`;
        } else {
            kVisit~=`
                if (newK.maxInnerDim==1){
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                } else {
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec1)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                }
                newK`~uniq~`=null;`;
        }
    } else if (pVisitors.length==9){
        if (pFlags==ParaFlags.Sequential){
            kVisit~=`
                if (!multiDim`~uniq~`){
                    if (newK.maxInnerDim==1){
                        newK`~uniq~`.exec0();
                    } else {
                        newK`~uniq~`.exec1();
                    }
                } else {
                    newK`~uniq~`.exec2();
                }`;
        } else {
            kVisit~=`
                if (!multiDim`~uniq~`){
                    if (newK.maxInnerDim==1){
                        Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec0)
                            .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                    } else {
                        Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec1)
                            .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                    }
                } else {
                    Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec2)
                        .appendOnFinish(&newK`~uniq~`.giveBack).autorelease.submitYield();
                }`;
        }
    } else {
        assert(0,"invalid number of elements in pVisitors");
    }
    kVisit~=`
            }
    `;
    res~=segArrayKLoopStr(pFlags,iterName,namesLocal,startLoop,kVisit,endLoop,uniq);
    res~="\n}\n";
    return res;
}

/// needs optimalBlockSize_1, and optimalBlockSize_2 (for loop 1, loop 2)
/// loop 1 (external)
/// public vars in kindVisit1: struct iterName, "context_1", "kIdx_1", "outOfRange_1", "visitKind_1"
/// public vars in pVisit1: name~"Ptr_1", "localPIndex_1", "pIndexPtr_1", "index_1"
/// loop 1 (internal, the vaiables in loop1 and)
/// public vars in kindVisit2: struct iterName, "context_2", "kIdx_2", "outOfRange_2", "visitKind_2"
/// public vars in pVisit2: name~"Ptr_2", "localPIndex_2", "pIndexPtr_2", "index_2"
char[] segArrayBinLoop(int pFlags,char[] iterName, char[][]namesLocal1, char[] contextExtra1,
    char[]startLoop1,char[]kindVisit1,char[]endLoop1,
    char[][]pVisitors1,
    char[][]namesLocal2, char[] contextExtra2,char[]startLoop2,char[] kindVisit2,
    char[]endLoop2,char[][]pVisitors2)
{
    char[] preLocalLoop2=`
        ParticleIdx index_1=outerLoopContext.lIndex;
        LocalPIndex localPIndex_1=LocalPIndex(outerLoopContext.kind,index_1);
        PIndex *pIndexPtr_1=outerLoopContext.pIndexPtrStart+cast(size_t)index_1;`;
    foreach(name;namesLocal1){
        preLocalLoop2~=`
        `~name~`.dtype *`~name~`Ptr_1=outerLoopContext.`~name~`PtrStart
                +outerLoopContext.`~name~`Nel*cast(size_t)outerLoopContext.lIndex;
        `;
    }
    char[][] pVis2;
    for(int i=0;i<pVisitors2.length;++i){
        if (i%3==0){
            pVis2~=preLocalLoop2~pVisitors2[i];
        } else {
            pVis2~=pVisitors2[i];
        }
    }
    auto intC=segArrayContextStr(pFlags,iterName~"_2", namesLocal2, contextExtra2~"\n"~iterName~"_1 *outerLoopContext;\n","_2",
    pVis2);
    
    char[][] pVis1;
    for(int i=0;i<pVisitors1.length;++i){
        if (i%3==1){
            pVis1~=`
            auto innerContext=innerLoopStartContext.alloc();
            innerLoopStartContext.outerLoopContext=this;
            `~pVisitors1[i]~`
            Task("particleLoop`~iterName~`",&innerContext.exec0).appendOnFinish(&innerContext.giveBack).autorelease.submitYield();
            `;
        } else {
            pVis1~=pVisitors1[i];
        }
    }

    auto extC=segArrayContextStr(pFlags,iterName~"_1", namesLocal1,contextExtra1~"\n"~intC~"\n"~iterName~"_2 *innerLoopStartContext;\n","_1",
        pVis1);
    auto kLoopInt=segArrayKLoopStr(pFlags,iterName~"_1."~iterName~"_2", namesLocal2,startLoop2,kindVisit2~`
    if (visitKind_2){
        newK_1.innerLoopStartContext=newK_2;
        newK_2.outerLoopContext=newK_1;
        auto tmpK1=mainContext_1.alloc();
        (*tmpK1)=(*newK_1);
        Task("kindIntLoop`~iterName~`",&newK_1.exec0).appendOnFinish(&newK_1.giveBack).autorelease.submitYield();
        newK_1=tmpK1;
    }
    `,endLoop2,"_2");
    auto kLoopExt=segArrayKLoopStr(pFlags,iterName~"_1", namesLocal1,startLoop1,kindVisit1~`
    if (visitKind_1){
        Task("kindMainLoop`~iterName~`",delegate void(){
        `~kLoopInt~`
        }).autorelease.submitYield();
    }
    `,endLoop1,"_1");
    
    char[] res=`{
    struct `~iterName~"_2;\n";
    res~=extC;
    res~=kLoopExt;
    res~="\n}\n";
    return res;
}
