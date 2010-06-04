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
        s(k.name)(" ");
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
    index_type[] kindStarts;
    T[] data;
    mixin(serializeSome("dchem.SegArray!("~T.stringof~")","kRange|kindStarts|data"));
    mixin printOut!();
    bool ownsData=true; // careful with positioning this, as it is the target of atomi ops...
    /// returns a writer that writes out the given array.
    /// memory is shared, so it reinteprets is as of type T
    static SegArrWriter opCall(V)(SegmentedArray!(V) sarr){
        SegArrWriter res;
        static assert((V.sizeof % T.sizeof)==0,"type of the array ("~V.stringof~") is not commensurate with "~T.stringof);
        if (sarr is null) return res;
        uint dimMult=V.sizeof/T.sizeof;
        if (dimMult!=1){
            res.kindStarts=new index_type[](sarr.kindStarts.length);
            foreach (i,p;sarr.kindStarts){
                res.kindStarts[i]=dimMult*p;
            }
        } else {
            res.kindStarts=sarr.kindStarts;
        }
        auto mdata=sarr.data.data;
        res.data=(cast(T*)mdata.ptr)[0..mdata.length*dimMult];
        res.kRange=sarr.kRange;
        res.ownsData=false;
        return res;
    }
    /// destroys the memory used by this object if it owns it
    void deallocData(){
        if (atomicCAS(ownsData,false,true)){
            delete data;
        }
    }
    /// returns a segmented array with this data as contents, if steal is true then the data will remain valid 
    /// also after the destruction of this. Memory is reinterpreted (no cast/conversion)
    SegmentedArray!(V) toSegArr(V=T)(SysStruct sysStruct,bool steal=false){
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
        return new SegmentedArray!(V)(aStruct,BulkArray!(V)((cast(V*)data.ptr)[0..data.length/dimMult],kRange),steal);
    }
    /// ditto
    SegmentedArray!(V) toSegArr(V=T)(SegmentedArrayStruct aStruct,bool steal=false,SegArrPool!(V) pool=null){
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
        auto newData=(cast(V*)data.ptr)[0..data.length/dimMult];
        BulkArray!(V) b;
        if (!atomicCAS(ownsData,false,true)) {
            b=BulkArray!(V)(newData.length);
            b.data()[]=newData;
        } else {
            b=BulkArray!(V)(newData);
        }
        return new SegmentedArray!(V)(aStruct,b,kRange,null,pool);
    }
    /// copies to an array of the given type, this might cast (checks kindStarts to decide if reinterpret the memory or cast/convert)
    void copyTo(V)(SegmentedArray!(V) s){
        assert(s!is null,"copy to only to allocated arrays");
        if (cast(void*)s.data.ptr is cast(void*)data.ptr) {
            assert(s.data.length*V.sizeof==data.length*T.sizeof,"unexpected length");
            return;
        }
        static if(is(T==V)){
            assert(kindStarts==s.kindStarts,"different kindStarts");
            s.data.data[]=data;
        } else static if (is(SegmentedArray!(V).basicDtype==T)){
            assert(data.length==s.data.basicData.length,"same basictype, but different lengths");
            s.data.basicData()[]=data;
        } else {
            if (kindStarts==s.kindStarts){
                selfArr=toSegArr(s.arrayStruct);
                s[]=selfArr;
            } else {
                if (s.data.length*V.sizeof==data.length*T.sizeof){
                    void*sPtr=s.data.ptr;
                    void*myPtr=data.ptr;
                    sPtr[0..data.length*T.sizeof]=myPtr[0..data.length*T.sizeof];
                } else {
                    assert(0,"could not copy to the array");
                }
            }
        }
    }
}

/// structure to dump out a DynPVector
struct DynPVectorWriter(T,int group){
    T[] cell;
    SegArrWriter!(T) pos;
    SegArrWriter!(T) orient;
    SegArrWriter!(T) dof;
    mixin(serializeSome("dchem.DynPVectorWriter!("~T.stringof~")","cell|pos|orient|dof"));
    mixin printOut!();
    /// true if this represents a non null the DynPVector
    bool isNonNull(){
        return pos.kindStarts.length!=0 || orient.kindStarts.length!=0 || dof.kindStarts.length!=0;
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
        if (v.cell!is null) res.cell=v.cell.h.cell;
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
            res.cell=new Cell(Matrix!(T,3,3)(cell),pVStruct.cellPeriod);
        }
        res.pos=pos.toSegArr!(Vector!(T,3))(pVStruct.posStruct,steal,pVStruct.poolPos);
        res.orient=orient.toSegArr!(Vector!(T,3))(pVStruct.orientStruct,steal,pVStruct.poolOrient);
        res.dof=dof.toSegArr!(Vector!(T,3))(pVStruct.dofStruct,steal,pVStruct.poolDof);
        return res;
    }
    /// copies to an array of the given type, this might cast (checks kindStarts to decide if reinterpret the memory or cast/convert)
    void copyTo(V,int group2)(ref DynPVectorWriter!(V,group2)v){
        static assert(group2==group,"copy only to same group");
        if (cell.length!=0){
            assert(cell.length==9,"invalid cell length");
            auto period=[0,0,0];
            if (v.cell is null){
                throw new Exception("cannot copy to null cell",__FILE__,__LINE__);
            }
            // create a new cell instead
            v.cell.h=Matrix!(T,3,3)(cell);
            v.cell.hInv=v.cell.h.inverse;
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
    DynPVectorWriter!(T,0) x;
    DynPVectorWriter!(T,1) dx;
    DynPVectorWriter!(T,1) mddx;
    Real potentialEnergy;
    Serializable hVars;
    mixin(serializeSome("dchem.PSysWriter!("~T.stringof~")","potentialEnergy|x|dx|mddx|hVars"));
    mixin printOut!();
    /// creates a writer for the given ParticleSys
    static PSysWriter opCall(ParticleSys!(T) pSys){
        PSysWriter res;
        res.potentialEnergy=pSys.potentialEnergy;
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
        pSys.potentialEnergy=potentialEnergy;
        // could be smarter about reusing memory...
        if (x.isNonNull()) {
            pSys.checkX();
            pSys.dynVars.x[]=x;
        }
        if (dx.isNonNull()){
            pSys.checkDx();
            pSys.dynVars.dx[]=dx;
        }
        if (mddx.isNonNull()){
            pSys.checkMddx();
            pSys.dynVars.mddx[]=mddx;
        }
        if (pSys.hVars !is null){
            pSys.hVars[]=hVars;
        }
    }
}

/// helper function, returns a struct that dumps the given ParticleSys
PSysWriter!(T) pSysWriter(T)(ParticleSys!(T) pSys){
    return PSysWriter!(T)(pSys);
}