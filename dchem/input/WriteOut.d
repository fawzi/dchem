module dchem.input.WriteOut;
import dchem.Physcon;
import blip.io.BasicIO;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.Common;
import dchem.sys.SegmentedArray;
import blip.serialization.Serialization;

/// writes out an xyz file
void writeXyz(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos,char[] comments){
    // angstrom*...
    auto nAtoms=sysStruct.particles[sysStruct.levels[0]].data.length;
    auto s=dumper(sink);
    s(nAtoms)("\n");
    foreach(c;comments){
        if (c=='\n' || c=='\r') throw new Exception("comments should have no newline",__FILE__,__LINE__);
    }
    s(comments)("\n");
    foreach(p;sysStruct.externalOrder.lSortedPIndex){
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        s(k.symbol)(" ");
        auto posAtt=pos[p,0];
        s(posAtt.x*angstrom)(" ")(posAtt.y*angstrom)(" ")(posAtt.z*angstrom)("\n");
    }
}

/// writes a turbomole style coordinate file
void writeTurboCoord(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos){
    auto s=dumper(sink);
    s("$coord\n");
    foreach(p;sysStruct.externalOrder.lSortedPIndex){
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        s(k.name)(" ");
        auto posAtt=pos[p,0];
        s(posAtt.x)(" ")(posAtt.y)(" ")(posAtt.z)(" ")(k.symbol)("\n");
    }
    s("$end\n");
}

/// structure to dump out a segmented array
struct SegArrWriter(T){
    KindRange kRange;
    uint[] kindStarts;
    T[] data;
    mixin(serializeSome("dchem.SegArray!("~T.stringof~")","kRange|kindStarts|data"));
    mixin printOut!();
    
    static SegArrWriter opCall(V)(SegmentedArray!(V) sarr){
        SegArrWriter res;
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        if (sarr is null) return res;
        uint dimMult=V.sizeof/T.sizeof;
        if (dimMult!=1){
            kindStarts=new index_type[](sarr.kindStarts.length);
            foreach (i,p;sarr.kindStarts){
                kindStarts[i]=dimMult*p;
            }
        } else {
            kindStarts=sarr.kindStarts;
        }
        auto mdata=sarr.data.data;
        data=(cast(T*)mdata.ptr)[0..mdata.length*dimMult];
        kRange=sarr.kRange;
    }
    
    SegmentedArray!(V) toSegArr(V=T)(SysStruct sysStruct){
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        uint dimMult=V.sizeof/T.sizeof;
        if (kindStarts.length==0) return null;
        kDims=new index_type[](kindStarts.length-1);
        auto pKStarts=sysStruct.particles.kindStarts;
        auto flags=SegmentedArrayStruct.Flags.Min1;
        foreach (i,p;kindStarts[1..$]){
            auto rFactor=(pKStarts[i+1]-pKStarts[i])*dimMult;
            assert(((p-kindStarts[i])%rFactor)==0,"number of elements non commensurate");
            kDims[i]=(p-kindStarts[i])/rFactor;
            if (kDims[i]==0) flags=0;
        }
        auto aStruct=new SegmentedArrayStruct(name,sysStruct.fullSystem,kRange,kDims,flags);
        return new SegmentedArray!(V)(aStruct,BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],kRange));
    }
    
    SegmentedArray!(V) toSegArr(V=T)(SegmentedArrayStruct aStruct){
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        uint dimMult=V.sizeof/T.sizeof;
        if (kindStarts.length==0) return null;
        auto pKStarts=sysStruct.particles.kindStarts;
        auto flags=SegmentedArrayStruct.Flags.Min1;
        foreach (i,p;kindStarts[1..$]){
            auto rFactor=(pKStarts[i+1]-pKStarts[i])*dimMult;
            if(((p-kindStarts[i])%rFactor)!=0)
                throw new Exception("number of elements non commensurate",__FILE__,__LINE__);
            if(aStruct.kindDim(kRange.kStart+i)!=(p-kindStarts[i])/rFactor)
                throw new Exception("different sizes",__FILE__,__LINE__);
            if (kDims[i]==0 && (aStruct.flags&SegmentedArrayStruct.Flags.Min1)!=0)
                throw new Exception("Min1 flags when not expected",__FILE__,__LINE__);
        }
        return new SegmentedArray!(V)(aStruct,BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],kRange));
    }
}

/// structure to dump out a DynPVector
struct DynPVectorWriter(T){
    T[] cell;
    SegArrWriter!(T) pos;
    SegArrWriter!(T) orient;
    SegArrWriter!(T) dof;
    mixin(serializeSome("dchem.DynPVectorWriter!("~T.stringof~")","cell|pos|orient|dof"));
    mixin printOut!();
    
    static DynPVectorWriter opCall(DynPVector!(T) v){
        DynPVectorWriter res;
        res.cell=v.cell.h.cell;
        res.pos=SegArrWriter!(T)(v.pos);
        res.orient=SegArrWriter!(T)(v.orient);
        res.dof=SegArrWriter!(T)(v.dof);
        return res;
    }
    
    DynPVector!(T) toDynPVector(V)(ParticleSys!(V) pSys){
        DynPVector!(T) res;
        if (cell.length!=0){
            assert(cell.length==9,"invalid cell length");
            auto period=[0,0,0];
            if (pSys.dynVars.x.cell !is null) period=pSys.dynVars.x.cell.period;
            res.cell=Cell(Matrix!(T,3,3)(cell),period);
        }
        if (pSys.dynVars.posStruct!is null){
            res.pos=pos.toSegArr!(Vector!(T,3))(pSys.dynVars.posStruct);
        } else {
            res.pos=pos.toSegArr!(Vector!(T,3))(pSys.sysStruct);
        }
        if (pSys.dynVars.orientStruct!is null){
            res.orient=orient.toSegArr!(Vector!(T,3))(pSys.dynVars.orientStruct);
        } else {
            res.orient=orient.toSegArr!(Vector!(T,3))(pSys.sysStruct);
        }
        if (pSys.dynVars.dofStruct!is null){
            res.dof=dof.toSegArr!(Vector!(T,3))(pSys.dynVars.dofStruct);
        } else {
            res.dof=dof.toSegArr!(Vector!(T,3))(pSys.sysStruct);
        }
    }
}
/// helper function to write out a DynPVector
DynPVectorWriter!(T) dynPVectorWriter(T)(DynPVector!(T) v){
    DynPVectorWriter!(T) res=DynPVectorWriter!(T)(v);
    return res;
}

/// structure to dump out a ParticleSystem
struct PSysWriter(T){
    DynPVectorWriter!(T) x;
    DynPVectorWriter!(T) dx;
    DynPVectorWriter!(T) mddx;
    Serializable hVars;
    mixin(serializeSome("dchem.PSysWriter!("~T.stringof~")","x|dx|mddx|hVars"));
    mixin printOut!();
    
    static PSysWriter opCall(ParticleSys!(T) pSys){
        PSysWriter res;
        res.x=dynPVectorWriter(pSys.x);
        res.dx=dynPVectorWriter(pSys.dx);
        res.mddx=dynPVectorWriter(pSys.mddx);
        res.hVars=pSys.hVars;
        return res;
    }
}
