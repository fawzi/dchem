module dchem.test.sys.SegmentedArrayTests;
import dchem.Common;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;
import blip.io.Console;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import dchem.test.sys.SubMappingTests;
import dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import dchem.sys.SubMapping;
import blip.sync.Atomic;

import blip.container.AtomicSLink;

struct RandomSegmentedArrayStruct(SegmentedArrayStruct.Flags flags=SegmentedArrayStruct.Flags.Min1,
    MappingKind mappingKind=MappingKind.Generic,char[]nameBase="RandomSArrayStruct"){
    SegmentedArrayStruct arrayStruct;
    
    static RandomSegmentedArrayStruct randomGenerate(Rand r, ref bool acceptable){
        RandomSegmentedArrayStruct res;
        RandomMap!(mappingKind,nameBase) submap;
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


struct RandomSegmentedArray(T,SegmentedArrayStruct.Flags flags=SegmentedArrayStruct.Flags.Min1,MappingKind mappingKind=MappingKind.Generic,char[]nameBase=""){
    SegmentedArrayStruct arrayStruct;
    Rand localR;
    SegmentedArray!(T)[] generatedArrays;
    
    static RandomSegmentedArray opCall(Rand r,SegmentedArrayStruct arrayStruct){
        RandomSegmentedArray res;
        res.localR=r;
        res.arrayStruct=arrayStruct;
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
        SegmentedArray!(T) res=new SegmentedArray!(T)(arrayStruct);
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
    sout("aStruct.arrayStruct.gSortedLocalPIndex.ptr:")(aStruct.arrayStruct.submapping.gSortedLocalPIndex.ptr)
        ("aStruct.arrayStruct.sortedPIndex.ptr:")(aStruct.arrayStruct.submapping.sortedPIndex.ptr)("\n");
    sout("loopTests start\n");
    auto sarr=new SegmentedArray!(LocalPIndex)(aStruct.arrayStruct);
    pragma(msg,segArrayMonoLoop(ParaFlags.Sequential,"iterContext",["array"],
        "int delegate(ref size_t i,ref PIndex, ref LocalPIndex, ref T) dlg; int* finalRes;","",
        "mainContext.dlg=loopBody;mainContext.finalRes=&result;","visitKind=visitKind&&(newK.pIndexPtr !is null);","",
        [`
        sout("localLoopIn\n");
        if ((*finalRes)!=0) return;
        sout("localLoopStart\n");
        `,"sout(\"localLoopIter{\\n\"); size_t ii=0; if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,*arrayPtr)){ *finalRes=returnV; return; } sout(\"localLoopIter}\\n\");","sout(\"localLoopEnd\\n\");",
        
        "sout(\"localLoop2In\\n\");if ((*finalRes)!=0) return;sout(\"localLoop2Start\\n\");",`
        sout("localLoop2Iter{\n");
        for (size_t ii=0;ii<arrayNel;++ii){
            if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
        }
        sout("localLoop2Iter}\n");
        `,"sout(\"localLoop2End\\n\");"]));
    
    sout("inline sLoop\n");
    {
        auto loopBody=delegate int(ref size_t i,ref PIndex pk, ref LocalPIndex lk,ref LocalPIndex v){
            sout("lpkIndexing\n");
            assert(sarr[lk,i]==v,"invalid indexing in sLoop2");
            sout("submap check\n");
            assert(sarr.arrayStruct.submapping[lk]==pk,"incorrect pk value");
            sout("pkIndexing\n");
            assert(sarr[pk][i]==v,"invalid indexing in sLoop3");
            sout("pkIndexingDone\n");
            return 0;
        };
        alias LocalPIndex T;
        auto array=sarr;
        size_t optimalBlockSize=s.val;
        int result;

        {

            struct iterContext{
                array.dtype* arrayPtrStart;
                size_t arrayNel;
                size_t optimalBlockSize;
                KindIdx kind;
                ParticleIdx start;
                ParticleIdx lIndex;
                ParticleIdx end;
                PIndex *pIndexPtrStart;
                size_t maxInnerDim;
                Exception exception;
                iterContext* context;
                iterContext* next;
        int delegate(ref size_t i,ref PIndex, ref LocalPIndex, ref T) dlg; int* finalRes;
                typeof(this) alloc(){
                    auto res=popFrom(context.next);
                    if (res is null){
                        res=new iterContext;
                    }
                    *res=*this;
                    res.next=this;
                    return res;
                }
                void exec0(){
                    if (context.exception !is null) return;

                sout("localLoopIn\n");
                if ((*finalRes)!=0) return;
                sout("localLoopStart\n");


                    {
                        try{
                            LocalPIndex localPIndex=LocalPIndex(kind,start);
                            PIndex * pIndexPtr=pIndexPtrStart;
                            array.dtype* arrayPtr=arrayPtrStart;
                            for (ParticleIdx index=start;index<end;++index){
                                lIndex=index;
                                sout("localLoopIter{\n"); size_t ii=0; if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,*arrayPtr)){ *finalRes=returnV; return; } sout("localLoopIter}\n");
                                ++localPIndex;
                                ++pIndexPtr;
                                arrayPtr+=arrayNel;
                            }
                        } catch (Exception e){
                            context.exception=e;
                        }
                    }
                    sout("localLoopEnd\n");
                }
                void exec1(){
                    if (context.exception !is null) return;
        sout("localLoop2In\n");if ((*finalRes)!=0) return;sout("localLoop2Start\n");

                    {
                        try{
                            LocalPIndex localPIndex=LocalPIndex(kind,start);
                            PIndex * pIndexPtr=pIndexPtrStart;
                            array.dtype* arrayPtr=arrayPtrStart;
                            for (ParticleIdx index=start;index<end;++index){
                                lIndex=index;

                sout("localLoop2Iter{\n");
                for (size_t ii=0;ii<arrayNel;++ii){
                    if (auto returnV=dlg(ii,*pIndexPtr,localPIndex,arrayPtr[ii])){ *finalRes=returnV; return; }
                }
                sout("localLoop2Iter}\n");

                                ++localPIndex;
                                ++pIndexPtr;
                                arrayPtr+=arrayNel;
                            }
                        } catch (Exception e){
                            context.exception=e;
                        }
                    }
                    sout("localLoop2End\n");
                }
                void giveBack(){
                    insertAt(context.next,this);
                }
            }

            {
                iterContext mainContext;
                auto kEnd=array.kRange.kEnd;
                iterContext* newK;
                PIndex dummyP;
                mainContext.optimalBlockSize=optimalBlockSize;
                mainContext.context=&mainContext;mainContext.dlg=loopBody;mainContext.finalRes=&result;
                    for (KindIdx kIdx=array.kRange.kStart;kIdx<kEnd;++kIdx){
                        if (mainContext.exception !is null) break;
                        bool visitKind=true;
                        bool outOfRange=false;
                        if (newK is null){
                            newK=mainContext.alloc();
                        }
                        newK.kind=kIdx;
                        newK.start=0;
                        newK.end=
                            array.arrayStruct.submapping.nLocalParticles(kIdx);
                        auto submap=array.arrayStruct.submapping;
                        if (submap.kindStarts[kIdx-submap.lKRange.kStart]<
                            submap.kindStarts[kIdx-submap.lKRange.kStart+1])
                        {
                            newK.pIndexPtrStart=submap.ptrI(
                                LocalPIndex(newK.kind,newK.start));
                        } else {
                            newK.pIndexPtrStart=null;
                        }
                        if (kIdx in array.kRange){
                            auto ik=cast(size_t)(kIdx-array.kRange.kStart);
                            if (array.kindStarts[ik]<array.kindStarts[ik+1]){
                                newK.arrayPtrStart=array._data.ptr+array.kindStarts[ik];
                                newK.arrayNel=array.arrayStruct.kindDim(kIdx);
                            } else {
                                outOfRange=true;
                                newK.arrayPtrStart=null;
                                newK.arrayNel=0;
                            }
                        } else {
                            outOfRange=true;
                            newK.arrayPtrStart=null;
                            newK.arrayNel=0;
                        }
        visitKind=visitKind&&(newK.pIndexPtrStart !is null);
                    if (visitKind){
                        newK.maxInnerDim=1;
                        bool multiDim=false;
                        if (newK.arrayNel!=newK.maxInnerDim && newK.arrayNel!=0 && newK.arrayNel!=1){
                            if (newK.maxInnerDim!=1){
                                multiDim=true;
                            }
                            if (newK.maxInnerDim<newK.arrayNel) newK.maxInnerDim=newK.arrayNel;
                        }
                        if (newK.maxInnerDim==1){
                            newK.exec0();
                        } else {
                            newK.exec1();
                        }
                        newK=null;
                    }


                }
                if (mainContext.exception !is null) throw new Exception("exception in SegmentedArray loop",__FILE__,__LINE__,mainContext.exception);
                auto freeL=mainContext.next;
                while (freeL!is null){
                    auto cNext=freeL.next;
                    delete freeL;
                    freeL=cNext;
                }

            }
        }
    }
    sout("calcRef\n");
    size_t iterRef=sarr.data.length;
    auto submap=sarr.arrayStruct.submapping;
    for (auto k=sarr.kRange.kStart;k<sarr.kRange.kEnd;++k){
        if (sarr.arrayStruct.kindDim(k)==0){
            iterRef-=1;
            iterRef+=submap.kindStarts[k-submap.lKRange.kStart+1]-submap.kindStarts[k-submap.lKRange.kStart];
        }
    }
    size_t ii=0;
    sout("sLoop\n");
    ii=0;
    foreach(i,ref v;sarr.sLoop){
        v=i;
        ++ii;
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");
    sout("sLoop2\n");
    foreach(i,v;sarr.sLoop){
        assert(v==i||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"invalid value in sLoop");
        assert(sarr[i,0]==v,"invalid indexing in sLoop");
    }
    sout("sLoop3\n");
    foreach(i,pk,lk,v;sarr.sLoop){
        assert(sarr[lk,i]==v,"invalid indexing in sLoop2");
        assert(sarr[pk][i]==v,"invalid indexing in sLoop3");
        assert(sarr.arrayStruct.submapping[lk]==pk,"incorrect pk value");
    }
    sout("sLoop4\n");
    foreach(pk,lk,v;sarr.sLoop){
        assert(lk==v||(sarr.arrayStruct.kindDim(lk.kind)==0 && v.kind==lk.kind),"value in sLoop");
        assert(sarr.arrayStruct.submapping[lk]==pk,"incorrect pk value");
    }
    sout("sLoop5\n");
    foreach(ref v;sarr.sLoop){
        auto lArr=sarr[v];
        assert(&v>=lArr.ptr && &v<lArr.ptr+lArr.length,"invalid indexing in sLoop");
    }
    sout("sLoop6\n");
    foreach(ref v;sarr.sLoop){
        v+=1;
    }
    sout("sLoop7\n");
    foreach(i,v;sarr.sLoop){
        assert(v.data==i.data+1||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"invalid value in sLoop");
    }
    
    sout("pLoop\n");
    ii=0;
    foreach(i,v;sarr.pLoop(s.val)){
        assert(v.data==i.data+1||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"incorrect value in pLoop");
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");
    sout("pLoop1\n");
    ii=0;
    foreach(i,pk,lk,v;sarr.pLoop(s.val)){
        assert(sarr[lk,i]==v,"invalid indexing in pLoop2");
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");
    sout("pLoop2\n");
    ii=0;
    foreach(pk,lk,v;sarr.pLoop(s.val)){
        assert(lk.data+1==v.data||(sarr.arrayStruct.kindDim(lk.kind)==0 && v.kind==lk.kind),"incorrect value in pLoop");
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");
    sout("pLoop3\n");
    ii=0;
    foreach(ref v;sarr.pLoop(s.val)){
        v+=1;
        atomicAdd!(size_t)(ii,1);
    }
    assert(ii==iterRef,"incorrect number of iterations in pLoop");
    sout("pLoop4\n");
    foreach(i,v;sarr.pLoop(s.val)){
        assert(i.data+2==v.data||(sarr.arrayStruct.kindDim(i.kind)==0 && v.kind==i.kind),"incorrect value in pLoop");
    }
    sout("loopTests end\n");
}


/// all container tests (a template to avoid compilation and instantiation unless really requested)
TestCollection segmentedArrayTests()(TestCollection superColl=null){
    TestCollection coll=new TestCollection("segmentedArray",__LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("testLoop",&testLoop,__LINE__,__FILE__,coll);
    return coll;
}

