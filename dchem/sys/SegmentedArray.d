module dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import blip.container.BulkArray;
import dchem.sys.SubMapping;
import dchem.Common;
import blip.narray.NArray;
import tango.core.sync.Atomic;
import tango.core.Traits;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import tango.core.Variant;

/// structure of a segmented array, kindStarts is valid only after freezing
final class SegmentedArrayStruct{
    enum Flags{
        None=0,   /// no special flags
        Frozen, /// kind dimensions cannot be changed anymore (this is set when an array is created)
        Direct, /// =submapping.mappingKind & MappingKind.Direct, worth caching?
        Min1,   /// at least one value per kind should be stored
    }
    SubMapping submapping;
    KindRange  kRange;
    index_type[]   _kindDims;
    index_type[]   kindStarts;
    Flags flags;
    mixin(serializeSome("dchem.sys.SegmentedArrayStruct","submapping|kRange|_kindDims|kindStarts"));
    // internal for serialization
    this(){ }
    /// allocates a new SegmentedArray with the given kind dimensions
    /// if the Min1 flag is set (default) then at least one value per kind is stored 
    this(SubMapping submapping,KindRange kRange,index_type[]kindDims,Flags f=Flags.Min1){
        assert(submapping !is null,"submapping needed to allocate");
        this.submapping = submapping;
        this.kRange     = kRange    ;
        this._kindDims   = kindDims  ;
        this.flags      = f;
        this.flags=(f& ~Flags.Direct)|(((submapping.mappingKind & MappingKind.Direct)!=0)?
                Flags.Direct:Flags.None);
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        assert(_kindDims.length==nkinds);
        assert(kRange in submapping.lKRange,"submapping is smaller than current kRange");
        if ((flags & Flags.Frozen)!=0){
            recalculateStarts();
        } else {
            kindStarts[]=index_type.max;
        }
    }
    /// ditto
    this(SubMapping submapping, KindRange kRange, NArray!(index_type,1) kindDims,Flags f=Flags.Min1){
        index_type[] kDims=new index_type[kindDims.shape[0]];
        for (size_t i=0;i<kindDims.shape[0];++i){
            kDims[i]=kindDims[i];
        }
        this(submapping,kRange,kDims,f);
    }
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
    index_type kindDim(KindIdx k){
        assert(k in kRange,"kind out of range");
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
        return new SegmentedArrayStruct(submapping,kRange,_kindDims.dup,flags & (~Flags.Frozen));
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
    
    mixin(serializeSome("dchem.sys.SegmentedArray","kRange|kindStarts|_data"));
    
    BulkArray!(T) data(){
        return _data[kindStarts[0],kindStarts[$]];
    }

    // internal for serialization
    this(){ }
    /// allocates a new SegmentedArray with the given kind dimensions
    /// min1 
    this(SegmentedArrayStruct arrayStruct, BulkArray!(T)data=BulkArray!(T).dummy,
        KindRange kRange=KindRange.all,index_type[] kindStarts=null)
    in{
        if (kindStarts!=null){
            auto nkinds=cast(size_t)(this.kRange.kEnd-this.kRange.kStart);
            auto kindShift=cast(size_t)(this.kRange.kStart-arrayStruct.kRange.kStart);
            for (size_t i=0;i<=nkinds;++i){
                assert(kindStarts[i]==kindStarts[0]+arrayStruct.kindStarts[i]-arrayStruct.kindStarts[kindShift],
                    "invalid kindStarts");
            }
        }
    } body {
        assert(arrayStruct is null,"arrayStruct must be valid");
        this.arrayStruct=arrayStruct.freeze;
        this.kRange=kRange;
        if (kRange.kEnd==KindIdx.init){
            assert(kRange.kStart==0);
            this.kRange=arrayStruct.kRange;
        }
        assert(this.kRange in arrayStruct.kRange);
        this.kindStarts=kindStarts;
        if (kindStarts==null){
            auto nkinds=cast(size_t)(this.kRange.kEnd-this.kRange.kStart);
            kindStarts=new index_type[cast(size_t)(nkinds+1)];
            auto kindShift=cast(size_t)(this.kRange.kStart-arrayStruct.kRange.kStart);
            for (size_t i=0;i<=nkinds;++i){
                kindStarts[i]=arrayStruct.kindStarts[kindShift+i]-arrayStruct.kindStarts[kindShift];
            }
        }
        _data=data;
        if (BulkArrayIsDummy(data)){
            _data=BulkArray!(T)(kindStarts[kindStarts.length-1]-kindStarts[0]+1);
        }
        direct=(arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0;
    }
    /// returns the array dimension of each particle of the given kind
    index_type kindDim(KindIdx kindIdx){
        return arrayStruct.kindDim(kindIdx);
    }
    /// array of elements in the given range
    SegmentedArray!(T) opIndex(KindRange kr){
        KindRange krCommon=kRange.intersect(kr);
        size_t startLK=cast(size_t)(krCommon.kStart-kRange.kStart);
        size_t endLK=cast(size_t)(krCommon.kEnd-kRange.kStart);
        return new SegmentedArray(arrayStruct,data,krCommon,kindStarts[startLK..endLK]);
    }
    /// array of elements for kind k
    BulkArray!(T) opIndex(KindIdx k){
        auto res=data[kindStarts[cast(size_t)k],kindStarts[cast(size_t)k+1]];
        // res.blockSize=kindDim[cast(size_t)k];
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
        auto startIdx=kindStarts[cast(size_t)(kindIdx-kRange.kStart)]+p.particle();
        return _data.getSlice(startIdx,startIdx+arrayStruct.kindDim(kindIdx));
    }
    /// array of elements for a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T[] opIndex(PIndex p){
        if ((arrayStruct.flags&SegmentedArrayStruct.Flags.Direct)!=0){
            auto kindIdx=p.kind();
            auto pos=p.particle();
            assert(kindIdx in kRange,"kind out of range");
            auto startIdx=kindStarts[cast(size_t)kindIdx]+p.particle();
            if (startIdx<kindStarts[cast(size_t)kindIdx+1])
                return _data.getSlice(startIdx,startIdx+arrayStruct.kindDim(kindIdx));
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            auto startIdx=kindStarts[cast(size_t)kindIdx]+l.particle();
            return _data.getSlice(startIdx,startIdx+arrayStruct.kindDim(kindIdx));
        }
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T *ptrI(LocalPIndex p,index_type i){
        auto kindIdx=p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange");
        assert(i>=0 && i<arrayStruct.kindDim(p.kind()),"index i out of bounds");
        auto startIdx=kindStarts[cast(size_t)(kindIdx-kRange.kStart)]+cast(size_t)p.particle();
        return _data.ptrI(startIdx+cast(size_t)i);
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    /// i is the index within the elements for particle i
    T opIndex(LocalPIndex p,index_type i){
        return *ptrI(p,i);
    }
    /// sets the value for element i of the particle using the local particle numbering
    /// (has to be in the kind range)
    void opIndexAssign(T val,LocalPIndex p,index_type i){
        *ptrI(p,i)=val;
    }
    /// address of element i of a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T *ptrI(PIndex p,index_type i){
        assert(i>=0 && i<arrayStruct.kindDim(p.kind()),"index i out of bounds");
        if (direct){
            auto kindIdx=p.kind();
            auto pos=p.particle();
            assert(kindIdx in kRange,"kind out of range");
            auto startIdx=kindStarts[cast(size_t)kindIdx]+p.particle();
            return _data.ptrI(startIdx+cast(size_t)i);
        } else {
            LocalPIndex l=arrayStruct.submapping[p];
            auto kindIdx=l.kind();
            assert(kindIdx in kRange,"kind out of range");
            auto startIdx=kindStarts[cast(size_t)kindIdx]+l.particle();
            return _data.ptrI(startIdx+cast(size_t)i);
        }
    }
    /// element i of a (global) particle p array
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T opIndex(PIndex p,index_type i){
        return *ptrI(p,i);
    }
    /// copies from an array to this
    void opSliceAssign(SegmentedArray val){
        assert(arrayStruct==val.arrayStruct,"different structs");
        assert(kRange==val.kRange,"different kRanges");
        assert(kindStarts==val.kindStarts,"different kindStarts");
        _data[]=val._data;
    }
    /// returns a copy of the segmented array
    SegmentedArray dup(){
        return new SegmentedArray(arrayStruct,_data.dup,kRange,kindStarts.dup);
    }
    /// returns a copy of the segmented array
    SegmentedArray deepdup(){
        return new SegmentedArray(arrayStruct,_data.deepdup,kRange,kindStarts.dup);
    }
}

char[] segArrayContextStr(char[] iterName, char[][]namesLocal,char[] contextExtra,
    char[]pVisitLocalStart,char[]pVisit,char[]pVisitLocalEnd,char[] uniq="")
{
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
            if (context is null)
                res.context=this;
            return res;
        }
        void exec(){`;
    visitorStr~="\n";
    visitorStr~=pVisitLocalStart;
    visitorStr~="\n";
    visitorStr~=`
            if (end-start>optimalBlockSize*3/2){
                auto mid=(start+end)/cast(ParticleIdx)2;
                auto firstHalf=alloc();
                firstHalf.end=mid;
                Task("SegArrLoop1`~iterName~`",&firstHalf.exec).autorelease.submitYield;
                auto secondHalf=alloc();
                secondHalf.start=mid;`;
    foreach (name;namesLocal){
        visitorStr~=`
                secondHalf.`~name~`PtrStart+=`~name~`Nel*mid;`;
    }
        visitorStr~=`
                Task("SegArrLoop2`~iterName~`",&secondHalf.exec).autorelease.submit;
            } else {
                LocalPIndex localPIndex`~uniq~`=LocalPIndex(kind,start);
                PIndex * pIndexPtr`~uniq~`=pIndexPtrStart;`;
        foreach (name;namesLocal){
            visitorStr~=`
                `~name~`.dtype* `~name~`Ptr`~uniq~`;`;
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
            }
            `~pVisitLocalEnd~`
        }
        void giveBack(){
            insertAt(context.next,this);
        }
        void exec2(){
            exec();
            giveBack();
        }
    }
    `;
    return visitorStr;
}


char[] segArrayKLoopStr(char[] iterName, char[][]namesLocal,
    char[]startLoop,char[] loopBody,char[]endLoop,char[] uniq="")
{
    char[] visitorStr=`
    {
        `~iterName~` mainContext`~uniq~`;
        auto kEnd`~uniq~`=`~namesLocal[0]~`.kRange.kEnd;
        `~iterName~`* newK`~uniq~`;
        mainContext`~uniq~`.optimalBlockSize=optimalBlockSize`~uniq~`;`;
    visitorStr~=startLoop;
    visitorStr~=`
        for (KindIdx kIdx`~uniq~`=`~namesLocal[0]~`.kRange.kStart;kIdx`~uniq~`<kEnd`~uniq~`;++kIdx`~uniq~`){
            bool visitKind`~uniq~`=true;
            bool outOfRange`~uniq~`=false;
            if (newK`~uniq~` is null){
                newK`~uniq~`=mainContext`~uniq~`.alloc();
            }
            newK`~uniq~`.kind=kIdx`~uniq~`;
            newK`~uniq~`.start=0;
            newK`~uniq~`.end=`~namesLocal[0]~`.arrayStruct.submapping.nLocalParticles(kIdx`~uniq~`);
            newK`~uniq~`.pIndexPtrStart=`~namesLocal[0]~`.arrayStruct.submapping.ptrI(LocalPIndex(kIdx`~uniq~`,cast(ParticleIdx)0));`;
    foreach (name;namesLocal){
        visitorStr~=`
            if (kIdx`~uniq~` in `~name~`.kRange){
                auto ik`~uniq~`=cast(size_t)(kIdx`~uniq~`-`~name~`.kRange.kStart);
                if (`~name~`.kindStarts[ik`~uniq~`]<`~name~`.kindStarts[ik`~uniq~`]){
                    newK`~uniq~`.`~name~`PtrStart=`~name~`._data.ptr+`~name~`.kindStarts[ik`~uniq~`];
                    newK`~uniq~`.`~name~`Nel=`~name~`.kindDim(kIdx`~uniq~`);
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
    visitorStr~=endLoop;
    visitorStr~=`
    }`;
    return visitorStr;
}

/// needs "optimalBlockSize"~uniq
/// public vars in kindVisit: struct iterName, "context"~uniq, "kIdx"~uniq, "outOfRange"~uniq, "visitKind"~uniq
/// public vars in pVisit: name~"Ptr"~uniq, "localPIndex"~uniq, "pIndexPtr"~uniq, "index"~uniq
char[] segArrayMonoLoop(char[] iterName, char[][]namesLocal,
    char[] contextExtra,char[]startLoop,char[]kindVisit,char[]endLoop,
    char[]pVisitLocalStart,char[]pVisit,char[]pVisitLocalEnd,char[] uniq="")
{
    char[] res="{\n";
    res~=segArrayContextStr(iterName, namesLocal, contextExtra,pVisitLocalStart,pVisit,pVisitLocalEnd,uniq);
    res~=segArrayKLoopStr(iterName, namesLocal,startLoop,kindVisit~uniq~`
    if (visitKind`~uniq~`){
        Task("kindMainLoop`~iterName~`",&newK`~uniq~`.exec2).autorelease.submitYield();
        newK`~uniq~`=null;
    }
    `,endLoop,uniq);
    res~="\n}\n";
    return res;
}

/// loop 1 (external) needs optimalBlockSize_1
/// public vars in kindVisit1: struct iterName, "context_1", "kIdx_1", "outOfRange_1", "visitKind_1"
/// public vars in pVisit1: name~"Ptr_1", "localPIndex_1", "pIndexPtr_1", "index_1"
/// loop 1 (internal, the vaiables in loop1 and) needs optimalBlockSize_2
/// public vars in kindVisit2: struct iterName, "context_2", "kIdx_2", "outOfRange_2", "visitKind_2"
/// public vars in pVisit2: name~"Ptr_2", "localPIndex_2", "pIndexPtr_2", "index_2"
char[] segArrayBinLoop(char[] iterName, char[][]namesLocal1, char[] contextExtra1,
    char[]startLoop1,char[]kindVisit1,char[]endLoop1,
    char[]pVisitLocalStart1,char[]pVisit1,char[]pVisitLocalEnd1,
    char[][]namesLocal2, char[] contextExtra2,char[]startLoop2,char[] kindVisit2,
    char[]endLoop2,char[]pVisitLocalStart2,char[]pVisit2,char[]pVisitLocalEnd2)
{
    char[] preLocalLoop2=`
        ParticleIdx index_1=outerLoopContext.lIndex;
        LocalPIndex localPIndex_1=LocalPIndex(outerLoopContext.kind,index_1);
        PIndex *pIndexPtr_1=outerLoopContext.pIndexPtrStart+cast(size_t)index_1;`;
    foreach(name;namesLocal1){
        preLocalLoop2~=`
        `~name~`.dtype *`~name~`Ptr_1=outerLoopContext.`~name~`PtrStart
                +outerLoopContext.`~name~`Nel*cast(size_t)outerLoopContext.lIndex;
        `~pVisitLocalStart2~"\n";
    }
    auto intC=segArrayContextStr(iterName~"_2", namesLocal2, contextExtra2~"\n"~iterName~"_1 *outerLoopContext;\n",
       preLocalLoop2,pVisit2,pVisitLocalEnd2,"_2");
    auto extC=segArrayContextStr(iterName~"_1", namesLocal1,contextExtra1~"\n"~intC~"\n"~iterName~"_2 *innerLoopStartContext;\n",
        pVisitLocalStart1,`
        auto innerContext=innerLoopStartContext.alloc();
        innerLoopStartContext.outerLoopContext=this;
        `~pVisit1~`
        Task("particleLoop`~iterName~`",&innerContext.exec2).autorelease.submitYield();
        `,pVisitLocalEnd1,"_1");
    auto kLoopInt=segArrayKLoopStr(iterName~"_1."~iterName~"_2", namesLocal2,startLoop2,kindVisit2~`
    if (visitKind_2){
        newK_1.innerLoopStartContext=newK_2;
        newK_2.outerLoopContext=newK_1;
        auto tmpK1=mainContext_1.alloc();
        (*tmpK1)=(*newK_1);
        Task("kindMainLoop`~iterName~`",&newK_1.exec2).autorelease.submitYield();
        newK_1=tmpK1;
    }
    `,endLoop2,"_2");
    auto kLoopExt=segArrayKLoopStr(iterName~"_1", namesLocal1,startLoop1,kindVisit1~`
    if (visitKind_1){
        `~kLoopInt~`
    }
    `,endLoop1,"_1");
    
    char[] res=`{
    struct `~iterName~"_2;\n";
    res~=extC;
    res~=kLoopExt;
    res~="\n}\n";
    return res;
}
