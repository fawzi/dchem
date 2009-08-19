module dchem.sys.SegmentedArray;
import dchem.sys.ParticleIndexes;
import blip.containers.BulkArray;
import dchem.sys.SubMapping;

/// segmented (level/kinds) array (with kind uniform dimension)
class SegmentedArray(T){
    SubMapping submapping;
    KindRange  kRange;
    size_t[]   kindDims;
    size_t[]   kindStarts;
    BulkArray!(T) _data;
    bool direct; // =submapping.mappingKind & MappingKind.Direct, worth caching?
    
    BulkArray!(T) data(){
        return _data[kindStarts[0],kindStarts[$]];
    }

    /// low level constructor
    this(SubMapping submapping,KindRange kRange,size_t[]kindDims,size_t[]kindStarts,BulkArray!(T) data){
        this.submapping = submapping;
        this.kRange     = kRange    ;
        this.kindDims   = kindDims  ;
        this.kindStarts = kindStarts;
        this._data      = data      ;
    }
    /// allocates a new SegmentedArray with the given kind dimensions
    /// min1 
    this(SubMapping submapping, KindRange kRange, IndexArray kindDims,bool min1=true){
        assert(submapping !is null,"submapping needed to allocate");
        this.kindDims=kindDims;
        this.min1=min1;
        this.kRange=kRange;
        this.direct=(submapping.mappingKind & MappingKind.Direct)!=0;
        this.submapping=submapping;
        auto nkinds=cast(size_t)(kRange.kEnd-kRange.kStart);
        assert(kindDims.length==nkinds);
        assert(kRange in submapping.lKRange,"submapping is smaller than current kRange");
        kindStarts=IndexArray(nkinds+1);
        
        kindStarts[0]=cast(size_t)0;
        auto kindShift=cast(size_t)(kRange.kStart-submapping.lKRange.kStart)
        for (size_t i=0;i<nkinds;++i){
            if (min1 && kindDims[i]==0){
                kindStarts[i+1]=kindStarts[i]+1; // always keep a "kind owned" value
            } else {
                kindStarts[i+1]=kindStarts[i]
                    +kindDims[i]*(submapping.kindStarts[i+1+kindShift]-submapping.kindStarts[i+kindShift]);
            }
        }
        data=BulkArray!(T)(kindStarts[kindStarts.length-1]);
    }
    /// number of elements of the kind k
    size_t kindDim(KindIdx k){
        return kindDims[cast(size_t)k]
    }
    /// array of elements in the given range
    SegmentedArray!(T) opIndex(KindRange kr){
        KindRange krCommon=kRange.intersect(kr);
        size_t startLK=cast(size_t)(krCommon.kStart-kRange.kStart);
        size_t endLK=cast(size_t)(krCommon.kEnd-kRange.kStart);
        return submapping,krCommon,kindDims[kindStarts[startLK..endLK],;
    }
    /// array of elements for kind k
    BulkArray!(T) opIndex(KindIdx k){
        auto res=data[kindStarts[cast(size_t)k],kindStarts[cast(size_t)k+1]];
        // res.blockSize=kindDims[cast(size_t)k];
        return res;
    }
    /// gets the particle using the local particle numbering, but that might be outside the kind range
    T[] getMaybeInRange(LocalPIndex p){
        auto kindIdx=cast(size_t)p.kind();
        if (kindIdx in kRange){
            auto startIdx=kindStarts[kindIdx]+p.position();
            return _data.getSlice(startIdx,startIdx+kindDims[kindIdx]);
        }
        return null;
    }
    /// gets the particle using the local particle numbering (has to be in the kind range)
    T[] opIndex(LocalPIndex p){
        auto kindIdx=cast(size_t)p.kind();
        assert(kindIdx in kRange,"kind is not in range as expected, you might want to use getMaybeInRange")
        auto startIdx=kindStarts[cast(size_t)(kindIdx-kRange.kStart)]+p.position();
        return _data.getSlice(startIdx,startIdx+kindDims[kindIdx]);
    }
    /// array of elements for a (global) particle p
    /// not as efficient as it could be, if you know that submapping is gapless for the
    /// kinds included cast to LocalPIndex and use either getMaybeInRange or opIndex
    T[] opIndex(PIndex p){
        if (direct){
            auto kindIdx=cast(size_t)p.kind();
            auto pos=p.position();
            if (kindIdx in kRange){
                auto startIdx=kindStarts[kindIdx]+p.position();
                if (startIdx<kindStarts[kindIdx+1])
                    return _data.getSlice(startIdx,startIdx+kindDims[kindIdx]);
            }
            return null;
        } else {
            LocalPIndex l=submapping[p];
            auto kindIdx=cast(size_t)l.kind();
            if (kindIdx in kRange){
                auto startIdx=kindStarts[kindIdx]+l.position();
                return _data.getSlice(startIdx,startIdx+kindDims[kindIdx]);
            }
        }
        return null;
    }
    
}

template visitorLocalStr(T,char[] nameV,char[][]namesLocal){
    const char[] visitorStr="",
    foreach (name;namesLocal){
        visitorStr~=`
        `~T.stringof~`* `~name~`Ptr0;
        `~T.stringof~`* `~name~`Ptr;
        `~T.stringof~`* `~name~`Inc=null;
        SegmentedArray!(`~T.stringof~`) `~name~`Arr;`;
    }
    // the following does not shortcut... avoid updating the rest if one is false?
    visitorStr~=`
    bool visitKind`~nameV~`(KindIdx k){
        bool res=true;`
    foreach (name;namesLocal){
        visitorStr~=`
        {
            auto resAtt=k in `~name~`Arr.kRange;
            if (resAtt){
                auto ik=cast(size_t)(k-`~name~`Arr.kRange.kStart);
                `~name~`Ptr0=`~name~`Arr._data.ptr+`~name~`Arr.kindStarts[ik];
                `~name~`Ptr=Ptr0;
                `~name~`inc=cast(`~T.stringof~`*)(cast(size_t)`~name~`Arr.kindDims[ik]*T.sizeof);
                if (`~name~`Arr.kindStarts[ik+1]>`~name~`Arr.kindStarts[ik]) return true;
                return false;
            } else {
                `~name~`Ptr0=null;
                `~name~`Ptr =null;
                `~name~`inc=null;
            }
        }`;
    }
    visitorStr~=`
        return res;
    }

    bool visitParticle`~nameV~`(LocalPIdx l){`;
    foreach (name; namesLocal){
        visitorStr~=`
        `~name~`Ptr=Ptr0+cast(size_t)l.position;`;
    }
    visitorStr~=`
        return true;
    }
    `;
}
