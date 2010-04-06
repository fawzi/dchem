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
import blip.io.BasicIO;
import blip.container.GrowableArray;
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
    char[] name;
    KindRange lKRange;
    index_type[] kindStarts;
    BulkArray!(PIndex) sortedPIndex;
    BulkArray!(LocalPIndex) gSortedLocalPIndex;
    BulkArray!(PIndex) lSortedPIndex;
    MappingKind mappingKind;

    mixin(serializeSome("dchem.sys.SubMapping","name|lKRange|kindStarts|sortedPIndex|gSortedLocalPIndex|lSortedPIndex"));
    mixin printOut!();

    ParticleIdx nLocalParticles(KindIdx k){
        assert(k in lKRange);
        auto ik=cast(size_t)(k-lKRange.kStart);
        return cast(ParticleIdx)(kindStarts[ik+1]-kindStarts[ik]);
    }
    /// lowlevel constructor
    static SubMapping opCall(char[] name,BulkArray!(PIndex) sortedPIndex,
        BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,index_type[] kindStarts, KindRange lKRange,MappingKind mappingKind)
    {
        return new SubMapping(name,sortedPIndex,gSortedLocalPIndex,lSortedPIndex,kindStarts,lKRange,mappingKind);
    }
    
    /// low level constructor, avoid its use, so that switching to struct would be easy
    this(char[] name,BulkArray!(PIndex) sortedPIndex,BulkArray!(LocalPIndex) gSortedLocalPIndex,
        BulkArray!(PIndex) lSortedPIndex,index_type[] kindStarts, KindRange lKRange,MappingKind mappingKind)
    {
        this.name=name;
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
        return lSortedPIndex.ptrI(cast(size_t)kindStarts[cast(size_t)(l.kind-lKRange.kStart)]
                             +cast(size_t)l.particle);
    }
    /// maps a local index to a global one
    /// faster mapping, does not work if sortedLocalPIndex has gaps
    PIndex opIndex(LocalPIndex l)
    in {
        assert(l.kind in lKRange,collectAppender(delegate void(CharSink s){
            dumperP(s)("indexing of a local particle with kind out of range:")
                (l)(" not in ")(lKRange);
        }));
        assert((cast(size_t)l.particle)+cast(size_t)kindStarts[l.kind-lKRange.kStart] <
            cast(size_t)kindStarts[l.kind-lKRange.kStart+1],collectAppender(delegate void(CharSink s){
                dumperP(s)("indexing of a local particle out of range:")
                    (l)(" when that kind has ")
                    (kindStarts[l.kind-lKRange.kStart+1]-kindStarts[l.kind-lKRange.kStart])
                    ("particles");
            }));
    } body {
        return lSortedPIndex[cast(size_t)kindStarts[cast(size_t)(l.kind-lKRange.kStart)]
                             +cast(size_t)l.particle];
    }
    /// performs various (expensive) tests on the submap consistency
    void testSubmap(){
        foreach(i,ref pk;sortedPIndex.pLoop){
            assert(this[this[pk]]==pk,"identity pk->lk->pk not valid");
            if(i!=0){
                auto prevPk=*(&pk-1);
                assert(prevPk<pk,collectAppender(delegate void(CharSink s){
                    dumperP(s)("non strictly ordered sortedPIndex [")(i-1)("]=")(prevPk)
                        (",[")(i)("]=")(pk)("\n");
                }));
            }
        }
        foreach(i,ref lk;gSortedLocalPIndex.pLoop){
            assert(this[this[lk]]==lk,"identity lk->pk->lk not valid");
            assert(this[lk]==sortedPIndex[i],collectAppender(delegate void(CharSink s){
                dumperP(s)("gSortedLocalPIndex ordering does not recover sortedPIndex at ")(i)(": ")(lk)
                    ("->")(this[lk])(" vs ")(sortedPIndex[i])("\n");
            }));
        }
        foreach(i,ref lk;lSortedPIndex.pLoop){
            if (i!=0){
                auto prevLk= *(&lk-1);
                assert(this[lk].data>this[prevLk].data,collectAppender(delegate void(CharSink s){
                    dumperP(s)("lSortedPIndex is not ordered [")(i-1)("]=")(prevLk)
                        (",[")(i)("]=")(lk)("\n");
                }));
            }
        }
        if (lSortedPIndex.length>0){
            assert(this[lSortedPIndex[0]].kind>=lKRange.kStart,"wrong lKRange (kStart too large)");
            assert(this[lSortedPIndex[lSortedPIndex.length-1]].kind<lKRange.kEnd,"wrong lKRange (kEnd too small)");
        }
    }
}
