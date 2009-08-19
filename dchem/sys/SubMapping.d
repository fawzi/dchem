module dchem.sys.SubMapping;
import dchem.sys.ParticleIndexes;
import blip.containers.BulkArray;

enum MappingKind:uint{
    Generic=0u,   /// generic mapping
    SameOrder=1u, /// order is manteined: sortedPIndex and lSortedPIndex are equal
    Gapless=2u,   /// whole ranges are mapped
    Direct=3u,    /// there is only a subsetting
    Same=7u,      /// no mapping at all
}
/// mapping os a system to a subsystem, assumes that the subsystem is smaller
/// (but only memory wastage occurs if this is not true)
class SubMapping{
    KindRange lKRange;
    size_t[] kindStarts;
    BulkArray!(PIndex) sortedPIndex;
    BulkArray!(LocalPIndex) gSortedLocalPIndex;
    BulkArray!(PIndex) lSortedPIndex;
    MappingKind mappingKind;

    /// lowlevel constructor
    static SubMapping opCall(BulkArray!(PIndex) sortedPIndex,
        BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,KindIdx kindStarts, KindRange lKRange){
    {
        return new SubMapping(sortedPIndex,gSortedLocalPIndex,lSortedPIndex,kindStarts,lKRange);
    }
    
    /// low level constructor, avoid its use, so that switching to struct would be easy
    this(BulkArray!(PIndex) sortedPIndex,BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,KindIdx kindStarts, KindRange lKRange){
        this.sortedPIndex=sortedPIndex;
        this.gSortedLocalPIndex=gSortedLocalPIndex;
        this.lSortedPIndex=lSortedPIndex;
        this.lKRange=lKRange;
        this.kindStarts=kindStarts;
        assert(kindStarts.length==1+cast(size_t)(lKRange.kEnd-lKRange.kStart));
        for (size_t i=1;i<kindStarts.length;++i){
            assert(kindStarts[i]<=kindStarts[i+1]);
        }
        assert(kindStarts.length>0);
        
    }
    
    /// maps a global index to a local one
    LocalPIndex opIndex(PIndex p){
        auto pos=sortedPIndex.data.lbound(p);
        if (sortedPIndex.length>pos && sortedPIndex[pos]==p){
            return gSortedLocalIndex[pos];
        }
        return LocalPIndex();
    }
    
    /// maps a local index to a global one
    /// faster mapping, does not work if sortedLocalPIndex has gaps
    PIndex opIndex(LocalPIndex l){
        assert(l in lKRange,"indexing of a local particle not in range");
        return lSortedPIndex[cast(size_t)kindStarts[cast(size_t)(l.kind-lKRange.kStart)]
                             +cast(size_t)l.position];
    }
    // evaluate in parallel
    class EvalOn{
        static Pool!(EvalOn) pool;
        static this(){
            pool=newPool!(EvalOn)
        }
        SubMapping s;
        LocalPIndex lPart;
        size_t lb,ub;
        Visitor v;
        this(s,lb,ub,v);
        EvalOn clear(){
            s=null;
            lPart=LocalPIndex();
            lb=size_t.max;
            ub=size_t.max;
            v.clear();
            return this;
        }
        void evaluate(){
            if (ub-lb>maxIdealBlockSize){
                auto mid=(ub+lb)/2;
                EvalOn(s,lb,mid,v).spawn()
                EvalOn(s,lb,mid,v).spawn()
            } else {
                LocalPIndex lP=lPart
                static if(is(typeof(visitor.visitParticle(size_t.init,LocalPIndex.init,PIndex.init)))){
                    for (size_t i=lb;i<ub;++i){
                        lP.data++;
                        visitor.visitParticle(i,lP,s.lSortedPIndex[i]);
                    }
                } else static if(is(typeof(visitor.visitParticle(size_t.init,LocalPIndex.init)))){
                    for (size_t i=lb;i<ub;++i){
                        lP.data++;
                        visitor.visitParticle(i,lP);
                    }
                } else {
                    static assert(false,"no low looping construct found");
                }
                static if(is(typeof(visitor.localEnd()))){
                    visitor.localEnd();
                }
            }
        }
    }
    // looping constructs
    template PVisitor(Visitor)(size_t maxIdealBlockSize=defaultBlockSize()){
        Visitor v;
        v.init()
        KindIdx k,kEnd;
        kEnd=lKRange.kEnd;
        for (k=lKRange.kStart;k<kEnd;++k){
            auto ik=k-lKRange.kStart;
            if (v.visitKind(k)){
                auto lb=kindStarts[ik];
                auto ub=kindStarts[ik+1];
                if (ub-lb>maxIdealBlockSize){
                    auto mid=(ub+lb)/2;
                    
                    visitor.spawn();
                }
            }
        }
    }
    
    
}
