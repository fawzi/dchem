module dchem.input.WriteOut;
import dchem.Physcon;
import blip.io.BasicIO;
import dchem.sys.DynVars;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.Common;
import dchem.sys.SegmentedArray;
import blip.serialization.Serialization;
import blip.sync.Atomic;
import blip.container.BulkArray;
import blip.util.Grow;
import stdlib=blip.stdc.stdlib;
import dchem.sys.Cell;

/// writes out an xyz file
void writeXyz(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos,char[] comments){
    // angstrom*...
    auto nAtoms=sysStruct.particles[sysStruct.levels[0]].length;
    auto s=dumper(sink);
    s(nAtoms)("\n");
    foreach(c;comments){
        if (c=='\n' || c=='\r') throw new Exception("comments should have no newline",__FILE__,__LINE__);
    }
    s(comments)("\n");
    foreach(p;sysStruct.externalOrder.gSortedLocalPIndex){
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        s(k.name)(" ");
        auto posAtt=pos[PIndex(p),0];
        s(posAtt.x*angstrom)(" ")(posAtt.y*angstrom)(" ")(posAtt.z*angstrom)("\n");
    }
}

/// writes a turbomole style coordinates.
/// skips the last newline
void writeTurboCoord(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos){
    auto s=dumper(sink);
    bool newLine=false;
    foreach(p;sysStruct.externalOrder.gSortedLocalPIndex){
        if (newLine) s("\n");
        newLine=true;
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        // s(k.name)(" ");
        auto posAtt=pos[PIndex(p),0];
        s(posAtt.x)(" ")(posAtt.y)(" ")(posAtt.z)(" ")(k.symbol);
    }
}

/// structure to dump out a segmented array
struct SegArrWriter(T){
    KindRange kRange;
    index_type[] kindStarts;
    index_type[] kindOffsets;
    BulkArray!(T) data;
    bool ownsData=true; // careful with positioning this, as it is the target of atomi ops...

    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(*this))("SegArrWriter!("~T.mangleof~")");
        metaI.addFieldOfType!(KindRange)("kRange","the range of this segment");
        metaI.addFieldOfType!(index_type[])("kindStarts","the starts of each kind");
        metaI.addFieldOfType!(LazyArray!(T))("data","the actual data"); // use T[]???
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    /// loops on the data (switch to pLoop?)...
    int opApply(int delegate(ref T el) loopBody){
        for (size_t ikind=0;ikind<kindOffsets.length;++ikind){
            auto sliceAtt=data[kindOffsets[ikind]..kindOffsets[ikind]+(kindStarts[ikind+1]-kindStarts[ikind])];
            auto resAtt=sliceAtt.opApply(loopBody);
            if (resAtt!=0) return resAtt;
        }
        return 0;
    }
    void serialize(Serializer s){
        s.field(metaI[0],kRange);
        s.field(metaI[1],kindStarts);
        auto dataS=LazyArray!(T)(&opApply,((kindStarts.length==0)?0:kindStarts[$-1]-kindStarts[0]));
        s.field(metaI[2],dataS);
    }
    void unserialize(Unserializer s){
        s.field(metaI[0],kRange);
        s.field(metaI[1],kindStarts);
        size_t len=0,startLen=8;
        if (kindStarts.length>0) startLen=kindStarts[$-1]-kindStarts[0];
        T[] dataA=(cast(T*)stdlib.malloc(startLen*T.sizeof))[0..startLen];
        auto dataS=LazyArray!(T)(delegate void(T el){
            if (len==startLen){
                startLen=growLength(startLen+1,T.sizeof);
                dataA=(cast(T*)stdlib.realloc(dataA.ptr,startLen*T.sizeof))[0..startLen];
            }
            dataA[len]=el;
            ++len;
        },delegate void(ulong l){
            if (l<startLen){
                startLen=l;
                dataA.length=startLen;
            }
        });
        s.field(metaI[2],dataS);
        auto guard=new ChunkGuard(dataA[0..len]); // leave it oversized??? then should make pool accept larger guards
        data=BulkArray!(T)(dataA[0..len],guard);
        kindOffsets=kindStarts;
        ownsData=true;
    }

    mixin printOut!();
    /// returns a writer that writes out the given array.
    /// memory is shared, so it reinteprets is as of type T
    static SegArrWriter opCall(V)(SegmentedArray!(V) sarr){
        SegArrWriter res;
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        if (sarr is null) return res;
        uint dimMult=V.sizeof/T.sizeof;
        auto supp=sarr.support;
        res.data=BulkArray!(T)((cast(T*)supp.ptr)[0..supp.length*dimMult],supp.guard);
        auto aStruct=sarr.arrayMap.arrayStruct;
        auto ikEnd=sarr.kRange.length;
        auto kShift=sarr.kRange.kStart-aStruct.kRange.kStart;
        if (dimMult!=1){
            res.kindStarts=new index_type[](sarr.kRange.length+1);
            for (size_t ik=0;ik<=ikEnd;++ik){
                res.kindStarts[ik]=aStruct.kindStarts[ik+kShift]*dimMult;
            }
        } else {
            res.kindStarts=aStruct.kindStarts[kShift..kShift+sarr.kRange.length+1];
        }
        res.kRange=sarr.kRange;
        res.kindOffsets=new index_type[](sarr.kRange.length);
        auto baseOffset=cast(size_t)(res.data.ptr-res.data.guard.dataPtr);
        for (size_t ik=0;ik<ikEnd;++ik){
            size_t nOffset=(sarr.kindOffsets[ik]-baseOffset)/T.sizeof;
            assert(((sarr.kindOffsets[ik]-baseOffset)%T.sizeof)==0,"non aligned segment");
            assert(nOffset+(res.kindStarts[ik+1]-res.kindStarts[ik])<=res.data.length,"offset overflow");
            res.kindOffsets[ik]=cast(index_type)nOffset;
        }
        res.kRange=sarr.kRange;
        res.ownsData=false;
        return res;
    }
    /// destroys the memory used by this object if it owns it
    void deallocData(){
        if (atomicCAS(ownsData,false,true)){
            if (data.guard!is null)
                data.guard.release();
        }
    }
    /// returns a segmented array with this data as contents, if steal is true then the data will remain valid 
    /// also after the destruction of this. Memory is reinterpreted (no cast/conversion)
    SegmentedArray!(V) toSegArrFromSysStruct(V=T)(SysStruct sysStruct,bool steal=false,char[] name="struct"){
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        uint dimMult=V.sizeof/T.sizeof;
        if (kindStarts.length==0) return null;
        auto kDims=new index_type[](kindStarts.length-1);
        auto pArray=sysStruct.particles;
        assert(kRange in pArray.kRange);
        auto pAStruct=pArray.arrayMap.arrayStruct;
        auto kOffset=kRange.kStart-pAStruct.kRange.kStart;
        auto pKStarts=pAStruct.kindStarts[kOffset..kOffset+kRange.kEnd+1];
        auto flags=SegmentedArrayStruct.Flags.Min1;
        foreach (i,p;kindStarts[1..$]){
            auto rFactor=(pKStarts[i+1]-pKStarts[i])*dimMult;
            assert(((p-kindStarts[i])%rFactor)==0,"number of elements non commensurate");
            kDims[i]=(p-kindStarts[i])/rFactor;
            if (kDims[i]==0) flags=SegmentedArrayStruct.Flags.None;
        }
        auto aStruct=new SegmentedArrayStruct(name,sysStruct.fullSystem,kRange,kDims,flags);
        auto aMap=new SegArrMemMap!(V)(aStruct);
        auto dData=BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],data.guard);
        if (steal && data.guard!is null){
            auto baseOffset=data.ptr-data.guard.dataPtr;
            bool compatible=true;
            auto kShift=kRange.kStart-aMap.kRange.kStart;
            foreach(i,o;kindOffsets){
                if (o*T.sizeof+baseOffset!=aMap.kindOffsets[i+kShift]) {
                    compatible=false;
                    break;
                }
            }
            // can be stolen:
            if (compatible && atomicCASB(ownsData,false,true)){
                return aMap.newArray(dData);
            }
        }
        auto res=aMap.newArray();
        foreach (k;kRange){
            auto ik=k-kRange.kStart;
            assert(kindOffsets[ik]%dimMult==0);
            assert((kindOffsets[ik]+kindStarts[ik+1]-kindStarts[ik])%dimMult==0);
            res[k][]=dData[kindOffsets[ik]/dimMult..
                (kindOffsets[ik]+kindStarts[ik+1]-kindStarts[ik])/dimMult];
        }
        return res;
    }
    /// ditto
    SegmentedArray!(V) toSegArr(V=T)(SegArrMemMap!(V) aMap,bool steal=false){
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        uint dimMult=V.sizeof/T.sizeof;
        if (kindStarts.length==0) return null;
        auto dData=BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],data.guard);
        // compatible layout:
        if (!(kRange in aMap.kRange))
            throw new Exception("incompatible kind ranges",__FILE__,__LINE__);
        if (steal && data.guard!is null){
            auto baseOffset=data.ptr-data.guard.dataPtr;
            bool compatible=true;
            auto kShift=kRange.kStart-aMap.kRange.kStart;
            foreach(i,o;kindOffsets){
                if (o*T.sizeof+baseOffset!=aMap.kindOffsets[i+kShift]) {
                    compatible=false;
                    break;
                }
            }
            // can be stolen:
            if (compatible && atomicCASB(ownsData,false,true)){
                return aMap.newArray(dData);
            }
        }
        auto res=aMap.newArray(kRange);
        foreach (k;kRange){
            auto ik=k-kRange.kStart;
            assert(kindOffsets[ik]%dimMult==0);
            assert((kindOffsets[ik]+kindStarts[ik+1]-kindStarts[ik])%dimMult==0);
            res[k][]=dData[kindOffsets[ik]/dimMult..
                (kindOffsets[ik]+kindStarts[ik+1]-kindStarts[ik])/dimMult];
        }
        return res;
    }
    /// copies to an array of the given type, this might cast (checks kindStarts to decide if reinterpret the memory or cast/convert)
    void copyTo(V)(SegmentedArray!(V) s){
        if (s is null) stdlib.abort();
        assert(s!is null,"copy to only to allocated arrays");
        uint dimMult=V.sizeof/T.sizeof;
        auto dData=BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],data.guard);
        foreach (k;kRange){
            auto ik=k-kRange.kStart;
            s[k][]=dData[kindOffsets[ik]*dimMult..
                (kindOffsets[ik]+kindStarts[ik+1]-kindStarts[ik])*dimMult];
        }
    }
}

/// structure to dump out a DynPVector
struct DynPVectorWriter(T,int group){
    T[] cell;
    T[] cellX0;
    SegArrWriter!(T) pos;
    SegArrWriter!(T) orient;
    SegArrWriter!(T) dof;
    int cellPeriod;
    mixin(serializeSome("dchem.DynPVectorWriter!("~T.stringof~")","cell|cellX0|cellPeriod|pos|orient|dof"));
    mixin printOut!();
    /// true if this represents a null DynPVector
    bool isDummy(){
        return !(pos.kindStarts.length!=0 || orient.kindStarts.length!=0 || dof.kindStarts.length!=0);
    }
    /// destroys the memory used by this object if it owns it
    void deallocData(){
        pos.deallocData();
        orient.deallocData();
        dof.deallocData();
    }
    /// greates a writer for the given vector
    static DynPVectorWriter opCall(DynPVector!(T,group) v){
        DynPVectorWriter res;
        if (v.cell!is null) {
            res.cell=v.cell.h.cell;
            res.cellX0=v.cell.x0.cell;
            res.cellPeriod=v.cell.periodicFlags;
        }
        res.pos=SegArrWriter!(T)(v.pos);
        res.orient=SegArrWriter!(T)(v.orient);
        res.dof=SegArrWriter!(T)(v.dof);
        return res;
    }
    /// returns a DynPVector of type V , if steal is true then the data will remain valid 
    /// also after the destruction of this. Memory is reinterpreted (no cast/conversion)
    DynPVector!(V,group) toDynPVector(V,int group2)(DynPVectorStruct!(V,group2) pVStruct,bool steal=false){
        static assert(group2==group,"conversion only within the same group");
        DynPVector!(T,group) res;
        if (cell.length!=0){
            assert(cell.length==9,"invalid cell length");
            Vector!(T,3) x0;
            if (cellX0.length!=0){
                assert(cellX0.length==3);
                x0=Vector!(T,3)(cellX0[0],cellX0[1],cellX0[2]);
            }
            res.cell=new Cell!(T)(Matrix!(T,3,3)(cell),cellPeriod,x0);
        }
        res.pos=pos.toSegArr!(Vector!(T,3))(pVStruct.poolPos,steal);
        res.orient=orient.toSegArr!(Vector!(T,3))(pVStruct.poolOrient,steal);
        res.dof=dof.toSegArr!(Vector!(T,3))(pVStruct.poolDof,steal);
        return res;
    }
    /// copies to an array of the given type, this might cast (checks kindStarts to decide if reinterpret the memory or cast/convert)
    void copyTo(V,int group2)(ref DynPVector!(V,group2)v){
        static assert(group2==group,"copy only to same group");
        if (cell.length!=0){
            assert(cell.length==9,"invalid cell length");
            Vector!(T,3) x0;
            if (cellX0.length!=0){
                assert(cellX0.length==3);
                x0=Vector!(T,3)(cellX0[0],cellX0[1],cellX0[2]);
            }
            v.cell=new Cell!(T)(Matrix!(T,3,3)(cell),cellPeriod,x0);
        }
        pos.copyTo(v.pos);
        orient.copyTo(v.orient);
        dof.copyTo(v.dof);
    }
}
/// helper function to write out a DynPVector
DynPVectorWriter!(T,group) dynPVectorWriter(T,int group)(DynPVector!(T,group) v){
    DynPVectorWriter!(T,group) res=DynPVectorWriter!(T,group)(v);
    return res;
}

/// structure to dump out a ParticleSystem
struct PSysWriter(T){
    DynPVectorWriter!(T,XType) x;
    DynPVectorWriter!(T,DxType) dx;
    DynPVectorWriter!(T,DxType) mddx;
    Real potentialEnergy;
    Real potentialEnergyError;
    Real mddxError;
    HiddenVars hVars;
    mixin(serializeSome("dchem.PSysWriter!("~T.stringof~")","potentialEnergy|potentialEnergyError|mddxError|x|dx|mddx|hVars"));
    mixin printOut!();
    /// creates a writer for the given ParticleSys
    static PSysWriter opCall(ParticleSys!(T) pSys){
        PSysWriter res;
        res.potentialEnergy=pSys.dynVars.potentialEnergy;
        res.potentialEnergyError=pSys.dynVars.potentialEnergyError;
        res.mddxError=pSys.dynVars.mddxError;
        res.x=dynPVectorWriter(pSys.dynVars.x);
        res.dx=dynPVectorWriter(pSys.dynVars.dx);
        res.mddx=dynPVectorWriter(pSys.dynVars.mddx);
        res.hVars=pSys.hVars;
        return res;
    }
    /// copies to an ParticleSys of the given type, this might cast 
    /// (checks kindStarts to decide if reinterpret the memory or cast/convert)
    void copyTo(V)(ParticleSys!(V) pSys){
        assert(pSys!is null,"cannot copy to null pSys");
        pSys.dynVars.potentialEnergy=potentialEnergy;
        pSys.dynVars.potentialEnergyError=potentialEnergyError;
        pSys.dynVars.mddxError=mddxError;
        // could be smarter about reusing memory...
        if (!x.isDummy()) {
            pSys.checkX();
            x.copyTo!(V,XType)(pSys.dynVars.x); // explicit instantiation of //pSys.dynVars.x[]=x;
        }
        if (!dx.isDummy()){
            pSys.checkDx();
            pSys.dynVars.dx[]=dx;
        }
        if (!mddx.isDummy()){
            pSys.checkMddx();
            pSys.dynVars.mddx[]=mddx;
        }
        if (pSys.hVars !is null){
            pSys.hVars.opSliceAssign(hVars);
        }
    }
    /// if the system is null
    bool isDummy(){
        return x.isDummy()&&dx.isDummy()&&mddx.isDummy()&&hVars is null;
    }
}

/// helper function, returns a struct that dumps the given ParticleSys
PSysWriter!(T) pSysWriter(T)(ParticleSys!(T) pSys){
    return PSysWriter!(T)(pSys);
}