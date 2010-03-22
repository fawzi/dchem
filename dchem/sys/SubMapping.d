/// submapping of a system, at the moment always explicit, but if Gapless
/// one could spare quite some memory, but then other code would need an if 
/// instead of always running on arrays
module dchem.sys.SubMapping;
import dchem.sys.PIndexes;
import blip.container.BulkArray;
import dchem.Common:index_type;
import tango.core.Array;
import blip.serialization.Serialization;
import blip.serialization.SerializationMixins;
import blip.BasicModels;
import blip.io.Console; // pippo

enum MappingKind:uint{
    Generic=0u,   /// generic mapping
    SameOrder=1u, /// order is manteined: sortedPIndex and lSortedPIndex are equal
    Gapless=2u,   /// whole ranges are mapped (if p1 and p2 are in then all p inbetween are in)
    Direct=3u,    /// there is only a subsetting
    Same=7u,      /// no mapping at all
}
/// mapping os a system to a subsystem, assumes that the subsystem is smaller
/// (but only memory wastage occurs if this is not true)
/// this object is supposed to be immutable once used...
class SubMapping: BasicObjectI{
    KindRange lKRange;
    index_type[] kindStarts;
    BulkArray!(PIndex) sortedPIndex;
    BulkArray!(LocalPIndex) gSortedLocalPIndex;
    BulkArray!(PIndex) lSortedPIndex;
    MappingKind mappingKind;

    mixin(serializeSome("dchem.sys.SubMapping","lKRange|kindStarts|sortedPIndex|gSortedLocalPIndex|lSortedPIndex"));
    mixin printOut!();

    ParticleIdx nLocalParticles(KindIdx k){
        assert(k in lKRange);
        auto ik=cast(size_t)(k-lKRange.kStart);
        return cast(ParticleIdx)(kindStarts[ik+1]-kindStarts[ik]);
    }
    /// lowlevel constructor
    static SubMapping opCall(BulkArray!(PIndex) sortedPIndex,
        BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,index_type[] kindStarts, KindRange lKRange,MappingKind mappingKind)
    {
        return new SubMapping(sortedPIndex,gSortedLocalPIndex,lSortedPIndex,kindStarts,lKRange,mappingKind);
    }
    
    /// low level constructor, avoid its use, so that switching to struct would be easy
    this(BulkArray!(PIndex) sortedPIndex,BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,index_type[] kindStarts, KindRange lKRange,MappingKind mappingKind)
    {
        this.sortedPIndex=sortedPIndex;
        this.gSortedLocalPIndex=gSortedLocalPIndex;
        this.lSortedPIndex=lSortedPIndex;
        this.lKRange=lKRange;
        this.kindStarts=kindStarts;
        this.mappingKind=mappingKind;
        sout("pippo kindStarts:")(kindStarts.length)(" vs ")(1+cast(size_t)(lKRange.kEnd-lKRange.kStart))("\n");
        assert(kindStarts.length==1+cast(size_t)(lKRange.kEnd-lKRange.kStart));
        for (size_t i=1;i<kindStarts.length;++i){
            assert(kindStarts[i-1]<=kindStarts[i]);
        }
        assert(kindStarts.length>0);
    }
    
    this(){ }
    
    /// maps a global index to a local one
    LocalPIndex opIndex(PIndex p){
        auto pos=lbound(sortedPIndex.data,p);
        if (sortedPIndex.length>pos && sortedPIndex[pos]==p){
            return gSortedLocalPIndex[pos];
        }
        return LocalPIndex.init;
    }
    
    /// returns a pointer into the lSortedPIndex
    PIndex *ptrI(LocalPIndex l){
        assert(l.kind in lKRange,"indexing of a local particle not in range");
        return &(lSortedPIndex[cast(size_t)kindStarts[cast(size_t)(l.kind-lKRange.kStart)]
                             +cast(size_t)l.particle]);
        
    }
    /// maps a local index to a global one
    /// faster mapping, does not work if sortedLocalPIndex has gaps
    PIndex opIndex(LocalPIndex l){
        assert(l.kind in lKRange,"indexing of a local particle not in range");
        return lSortedPIndex[cast(size_t)kindStarts[cast(size_t)(l.kind-lKRange.kStart)]
                             +cast(size_t)l.particle];
    }
    
    /+
    // evaluate in parallel
    class EvalOn(Visitor){
        Pool!(EvalOn) pool;
        SubMapping s;
        LocalPIndex lPart;
        size_t lb,ub;
        Visitor v;
        this(){}
        this(SubMapping s,LocalPIndex lb,LocalPIndex ub,Visitor v,Pool!(EvalOn) pool=null){
            this.s=s;
            this.lb=lb;
            this.ub=ub;
            this.v=v;
            this.pool=pool;
        }
        EvalOn clear(){
            s=null;
            lPart=LocalPIndex();
            lb=size_t.max;
            ub=size_t.max;
            v.clear();
            return this;
        }
        reset(SubMapping s,LocalPIndex lb,LocalPIndex ub,Visitor v,Pool!(EvalOn) pool=null){
            this.s=s;
            this.lb=lb;
            this.ub=ub;
            this.v=v;
            this.pool=pool;
        }
        void evaluate(){
            if (ub-lb>maxIdealBlockSize){
                auto mid=(ub+lb)/2;
                if (pool is null){
                    pool=new Pool!(EvalOn)();
                }
                Task(&(pool.getObj().reset(s,lb,mid,v,pool).evaluate)).autorelease.spawnYield();
                Task(&(pool.getObj().reset(s,lb,mid,v,pool).evaluate)).autorelease.spawn();
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
            if (this.pool!is null){
                pool.giveBack(this);
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
                    auto eval=EvalOn!(Visitor)(this,lb,ub,v);
                    eval.evaluate()
                } else {
                    LocalPIndex lP=LocalPIndex(ik,lb);
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
    }
    +/
}
