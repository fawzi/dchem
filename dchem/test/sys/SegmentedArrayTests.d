module dchem.test.sys.SegmentedArrayTests;
import dchem.Common;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import dchem.test.sys.SubMappingTests;
import dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import dchem.sys.SubMapping;
import blip.sync.Atomic;

struct RandomSegmentedArrayStruct(SegmentedArrayStruct.Flags flags=SegmentedArrayStruct.Flags.Min1,
    MappingKind mappingKind=MappingKind.KindPreserving,char[]nameBase="RandomSArrayStruct"){
    SegmentedArrayStruct arrayStruct;
    
    static RandomSegmentedArrayStruct randomGenerate(Rand r, ref bool acceptable){
        RandomSegmentedArrayStruct!(flags,mappingKind|MappingKind.KindPreserving) res;
        RandomMap!(mappingKind|MappingKind.KindPreserving,nameBase) submap;
        acceptable=acceptable&&simpleRandom(r,submap);
        char[] namePostfix;
        simpleRandom(r,namePostfix);
        KindRange  kRange;
        if (r.uniform!(bool)()){
            kRange.kStart=submap.map.lKRange.kStart;
        } else{
            kRange.kStart=cast(KindIdx)r.uniformR2(cast(int)submap.map.lKRange.kStart,
                cast(int)submap.map.lKRange.kEnd+1);
        }
        if (r.uniform!(bool)()){
            kRange.kEnd=submap.map.lKRange.kEnd;
        } else{
            kRange.kEnd=cast(KindIdx)(cast(int)kRange.kStart+
                r.uniformR(submap.map.lKRange.kEnd-kRange.kStart+1));
        }
        index_type[]   kindDims=new index_type[](kRange.kEnd-kRange.kStart);
        if (r.uniform!(bool)){
            kindDims[]=1;
        } else if (r.uniform!(bool)){
            foreach (ref kD; kindDims){
                kD=r.uniformR(2);
            }
        } else {
            foreach (kD; kindDims){
                kD=generateSize(r);
            }
        }
        res.arrayStruct=new SegmentedArrayStruct(nameBase~namePostfix,submap.map,kRange,kindDims,flags);
        return res;
    }
    void desc(CharSink c){
        writeOut(c,arrayStruct);
    }
}


struct RandomSegmentedArray(T,SegmentedArrayStruct.Flags flags=SegmentedArrayStruct.Flags.Min1,MappingKind mappingKind=MappingKind.KindPreserving,char[]nameBase=""){
    SegmentedArrayStruct arrayStruct;
    SegArrMemMap!(T) arrayMap;
    Rand localR;
    SegmentedArray!(T)[] generatedArrays;
    
    static RandomSegmentedArray opCall(Rand r,SegmentedArrayStruct arrayStruct){
        RandomSegmentedArray res;
        res.localR=r;
        res.arrayStruct=arrayStruct;
        res.arrayMap=new SegArrMemMap!(T)(arrayStruct);
        return res;
    }
    static RandomSegmentedArray randomGenerate(Rand r, ref bool acceptable){
        RandomSegmentedArray res;
        RandomSegmentedArrayStruct!(flags,mappingKind) randomAStruct;
        acceptable=acceptable && simpleRandom(r,randomAStruct);
        res.localR=r.spawn();
        res.arrayStruct=randomAStruct.arrayStruct;
        return res;
    }
    /// returns a random segmented array, null if it can't generate it
    SegmentedArray!(T) randomArray(){
        SegmentedArray!(T) res=arrayMap.newArray();
        for (int i=0;i<10;++i){
            bool acceptable=true;
            mkRandomArray(localR,res.data.data,acceptable);
            if (acceptable) {
                generatedArrays~=res;
                return res;
            }
        }
        return null;
    }
    void desc(CharSink c){
        c("{ arrayStruct:");
        writeOut(c,arrayStruct);
        c(" generatedArrays:");
        writeOut(c,generatedArrays);
        c("}");
    }
}

void testLoop(RandomSegmentedArrayStruct!() aStruct,SizeLikeNumber!(3,1) s){
    auto aMap=new SegArrMemMap!(LocalPIndex)(aStruct.arrayStruct);
    auto sarr=aMap.newArray();
    
    size_t iterRef=sarr.length;
    size_t ii=0;
    ii=0;
    foreach(i, vArr;sarr.sLoop){
        vArr[]=i;
        ++ii;
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");

    foreach(i,vArr;sarr.sLoop){
        foreach (ref v;vArr){
            assert(v==i||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"invalid value in sLoop");
            assert(sarr[i,0]==v,"invalid indexing in sLoop");
        }
    }

    foreach(i,pk,lk,v;sarr.sLoop){
        assert(sarr[lk,i]==v,"invalid indexing in sLoop2");
        assert(sarr[pk][i]==v,"invalid indexing in sLoop3");
        assert(sarr.arrayStruct.submapping[lk]==pk,"incorrect pk value");
    }

    foreach(pk,lk,vArr;sarr.sLoop){
        foreach(v;vArr){
            assert(lk==v||(sarr.arrayStruct.kindDim(lk.kind)==0 && v.kind==lk.kind),"value in sLoop");
            assert(sarr.arrayStruct.submapping[lk]==pk,"incorrect pk value");
        }
    }

    foreach(ref v;sarr.sLoop){
        auto lArr=sarr[v];
        assert(&v>=lArr.ptr && &v<lArr.ptr+lArr.length,"invalid indexing in sLoop");
    }

    foreach(ref v;sarr.sLoop){
        v+=1;
    }

    foreach(i,vArr;sarr.sLoop){
        foreach(v;vArr)
            assert(v.data==i.data+1||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"invalid value in sLoop");
    }
    
    foreach(i,vArr;sarr.sLoop){
        foreach(v;vArr){
            sarr[i,0]=LocalPIndex(sarr[i,0].data-1);
            assert(sarr[i,0]==LocalPIndex(v.data-1),"element set failed");
            sarr[i,0]=v;
        }
    }
    
    ii=0;
    foreach(i,vArr;sarr.pLoop(s.val)){
        foreach(v;vArr){
            assert(v.data==i.data+1||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"incorrect value in pLoop");
            atomicAdd!(size_t)(ii,1);
        }
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");

    ii=0;
    foreach(i,pk,lk,v;sarr.pLoop(s.val)){
        assert(sarr[lk,i]==v,"invalid indexing in pLoop2");
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");

    ii=0;
    foreach(pk,lk,vArr;sarr.pLoop(s.val)){
        foreach(v;vArr){
            assert(lk.data+1==v.data||(sarr.arrayStruct.kindDim(lk.kind)==0 && v.kind==lk.kind),"incorrect value in pLoop");
            atomicAdd!(size_t)(ii,1);
        }
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");

    ii=0;
    foreach(ref v;sarr.pLoop(s.val)){
        v+=1;
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");

    foreach(i,vArr;sarr.pLoop(s.val)){
        foreach(v;vArr)
            assert(i.data+2==v.data||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"incorrect value in pLoop");
    }
}

/// all container tests (a template to avoid compilation and instantiation unless really requested)
TestCollection segmentedArrayTests()(TestCollection superColl=null){
    TestCollection coll=new TestCollection("segmentedArray",__LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("testLoop",&testLoop,__LINE__,__FILE__,coll);
    return coll;
}

