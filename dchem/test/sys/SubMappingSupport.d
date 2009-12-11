/// generation of random submappings
module dchem.test.sys.SubMappingSupport;
import dchem.sys.SubMapping;
import tango.math.random.Random;

/// returns a submap that has at least the requested mappingKind
SubMapping randomSubmapping(Random r,SubMapping mainSubmap,MappingKind mappingKind){
    MappingKind[5] mapK=[MappingKind.Same,MappingKind.Direct,MappingKind.Gapless,
        MappingKind.SameOrder,MappingKind.Generic];
    MappingKind toGen;
    switch(mappingKind){
    case MappingKind.Generic:
        toGen=cast(MappingKind)r.uniformR1(5);
        break;
    case MappingKind.SameOrder:
        toGen=cast(MappingKind)r.uniformR1(3);
        if (toGen==MappingKind.Gapless) toGen=MappingKind.SameOrder;
        break;
    case MappingKind.Gapless:
        toGen=cast(MappingKind)r.uniformR1(3);
        break;
    case MappingKind.Direct:
        toGen=cast(MappingKind)r.uniformR1(2);
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
            auto nKindsMain=cast(int)(mainSubmap.lKRange.kEnd-mainSubmap.lKRange.kStart);
            auto nKinds=r.uniformR1(nKindsMain);
            res.lKRange.kStart=cast(KindIdx)(r.uniformR1(nKindsMain-nKindsMain+1);
            res.lKRange.kEnd=res.lKRange.kStart+cast(KindIdx)nKinds;
            res.kindStarts=new index_type[](nKinds+1);
            res.kindStarts[0]=0;
            for (int i=0;i<nKinds;++i){
                res.kindStarts[i+1]=res.kindStarts[i]
                    +r.uniformR1(mainSubmap.kindStarts[i+1]-mainSubmap.kindStarts[i]);
            }
            auto sortedPIndex=BulkArray!(PIndex)(kindStarts[$]);
            auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(kindStarts[$]);
            for (int i=0;i<nKinds;++i){
                auto nP=res.kindStarts[i+1]-res.kindStarts[i];
                rStart=r.uniformR1(mainSubmap.kindStarts[i+1]-mainSubmap.kindStarts[i]-nP);
                auto kind=res.lKRange.kStart+cast(KindIdx)i;
                auto pidx=PIndex(kind,cast(ParticleIdx)rStart);
                auto lidx=LocalPIndex(kind,cast(ParticleIdx)0);
                for (size_t p=0;p<nP;++p){
                    sortedPIndex[p]=pidx;
                    pidx+=1;
                }
                for (size_t p=0;p<nP;++p){
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
            return mainSubmap.dup;
        default:
            throw new Exception("Unkown mapping Kind requested",__FILE__,__LINE__);
    }
}
