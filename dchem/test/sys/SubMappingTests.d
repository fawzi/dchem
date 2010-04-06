/// generation of random submappings
module dchem.test.sys.SubMappingTests;
import dchem.sys.SubMapping;
import tango.math.random.Random;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;
import dchem.sys.PIndexes;
import blip.container.BulkArray;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.io.Console;

/// an identity mapping on the kinds 0..nParticlesPerKind.length, with each nParticlesPerKind[i]
/// particles
SubMapping identitySubmapping(char[] name,uint[] nParticlesPerKind){
    auto res=new SubMapping();
    res.name=name;
    auto nKinds=nParticlesPerKind.length;
    res.lKRange.kStart=cast(KindIdx)0;
    res.lKRange.kEnd=cast(KindIdx)nKinds;
    res.kindStarts=new index_type[](nKinds+1);
    res.kindStarts[0]=0;
    for (int i=0;i<nKinds;++i){
        res.kindStarts[i+1]=res.kindStarts[i]+nParticlesPerKind[i];
    }
    auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(res.kindStarts[$-1]);
    for (int i=0;i<nKinds;++i){
        auto kind=cast(KindIdx)i;
        auto lidx=LocalPIndex(kind,cast(ParticleIdx)0);
        for (size_t p=res.kindStarts[i];p<res.kindStarts[i+1];++p){
            gSortedLocalPIndex[p]=lidx;
            lidx+=1;
        }
    }
    auto sortedPIndex=BulkArray!(PIndex)((cast(PIndex*)gSortedLocalPIndex.ptr)[0..gSortedLocalPIndex.length],
        gSortedLocalPIndex.guard);
    res.sortedPIndex=sortedPIndex;
    res.gSortedLocalPIndex=gSortedLocalPIndex;
    res.lSortedPIndex=sortedPIndex;
    res.mappingKind=MappingKind.Same;
    return res;
}

/// returns a submap that has at least the requested mappingKind
/// to do: add a real generic mapping...
SubMapping randomSubmapping(char[] baseName,Random r,SubMapping mainSubmap,MappingKind mappingKind){
    MappingKind[5] mapK=[MappingKind.Same,MappingKind.Direct,MappingKind.Gapless,
        MappingKind.SameOrder,MappingKind.Generic];
    MappingKind toGen;
    char[] namePostfix;
    simpleRandom(r,namePostfix);
    auto name=baseName~namePostfix;
    switch(mappingKind){
    case MappingKind.Generic:
        toGen=mapK[r.uniformR(5)];
        break;
    case MappingKind.SameOrder:
        toGen=mapK[r.uniformR(3)];
        if (toGen==MappingKind.Gapless) toGen=MappingKind.SameOrder;
        break;
    case MappingKind.Gapless:
        toGen=mapK[r.uniformR(3)];
        break;
    case MappingKind.Direct:
        toGen=mapK[r.uniformR(2)];
        break;
    case MappingKind.Same:
        toGen=MappingKind.Same;
        break;
    default:
        throw new Exception("Unkown mapping Kind requested",__FILE__,__LINE__);
        break;
    }
    switch(toGen){
        case MappingKind.Generic:
        case MappingKind.SameOrder:
        case MappingKind.Gapless:
        case MappingKind.Direct:
            auto res=new SubMapping();
            res.name=name;
            auto nKindsMain=cast(int)(mainSubmap.lKRange.kEnd-mainSubmap.lKRange.kStart);
            auto nKinds=r.uniformR(nKindsMain+1);
            res.lKRange.kStart=cast(KindIdx)(r.uniformR(nKindsMain-nKinds+1));
            res.lKRange.kEnd=cast(KindIdx)(res.lKRange.kStart+nKinds);
            res.kindStarts=new index_type[](nKinds+1);
            res.kindStarts[0]=0;
            for (int i=0;i<nKinds;++i){
                res.kindStarts[i+1]=res.kindStarts[i]
                    +r.uniformR(mainSubmap.kindStarts[i+1]-mainSubmap.kindStarts[i]+1);
            }
            auto sortedPIndex=BulkArray!(PIndex)(res.kindStarts[$-1]);
            auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(res.kindStarts[$-1]);
            for (int i=0;i<nKinds;++i){
                auto nP=res.kindStarts[i+1]-res.kindStarts[i];
                auto rStart=r.uniformR(mainSubmap.kindStarts[i+1]-mainSubmap.kindStarts[i]-nP+1);
                auto kind=res.lKRange.kStart+cast(KindIdx)i;
                auto pidx=PIndex(kind,cast(ParticleIdx)rStart);
                auto lidx=LocalPIndex(kind,cast(ParticleIdx)0);
                for (size_t p=res.kindStarts[i];p<res.kindStarts[i+1];++p){
                    sortedPIndex[p]=pidx;
                    pidx+=1;
                }
                for (size_t p=res.kindStarts[i];p<res.kindStarts[i+1];++p){
                    gSortedLocalPIndex[p]=lidx;
                    lidx+=1;
                }
            }
            res.sortedPIndex=sortedPIndex;
            res.gSortedLocalPIndex=gSortedLocalPIndex;
            res.lSortedPIndex=sortedPIndex;
            if (r.uniform!(bool)()){
                res.mappingKind=mappingKind;
            } else {
                res.mappingKind=toGen;
            }
            return res;
        case MappingKind.Same:
            return mainSubmap;
        default:
            sout("kind:")(toGen);
            throw new Exception("Unkown mapping Kind requested",__FILE__,__LINE__);
    }
}

struct RandomMap(MappingKind mappingKind=MappingKind.Generic,char[] baseName=""){
    SubMapping map;
    static RandomMap randomGenerate(Rand r){
        RandomMap res;
        auto nKinds=generateSize(r);
        auto nParticlesPerKind=new uint[](nKinds);
        foreach(ref np;nParticlesPerKind){
            np=generateSize(r);
        }
        char[] namePostfix;
        simpleRandom(r,namePostfix);
        auto name=baseName~namePostfix;
        auto mainMap=identitySubmapping(name,nParticlesPerKind);
        switch(mappingKind){
            case MappingKind.Same:
            res.map=mainMap;
            break;
            case MappingKind.Generic:
            case MappingKind.SameOrder:
            case MappingKind.Gapless:
            case MappingKind.Direct:
            res.map=mainMap; // pippo to do make real random...
            //res.map=randomSubmapping(name,r,mainMap,mappingKind);
            break;
            default:
                throw new Exception("Unkown mapping Kind requested",__FILE__,__LINE__);
        }
        return res;
    }
    void desc(CharSink s){
        writeOut(s,map);
    }
}

void testRandomMaps(RandomMap!() rMap){
    rMap.map.testSubmap();
}

/// collection of all the tests on the submap module
TestCollection submapTests(TestCollection superColl){
    TestCollection coll=new TestCollection("submap",
        __LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("randomMapTest",&testRandomMaps,__LINE__,__FILE__,coll);
    return coll;
}

