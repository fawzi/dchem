/// converts from ReadIn to particle system and viceversa
module dchem.input.ReadIn2PSys;
import dchem.input.ReadIn;
import dchem.sys.ParticleSys;
import blip.util.NotificationCenter;
import dchem.sys.PIndexes;
import dchem.PeriodicTable;
import blip.container.BulkArray;
import blip.container.GrowableArray;
import blip.io.BasicIO;
import dchem.Common;
import dchem.sys.SubMapping;
import dchem.sys.SegmentedArray;
import blip.sync.Atomic;
import dchem.sys.Cell;
import blip.io.Console; // pippo

class ParticleKindMap{
    ParticleKind[char[]] specialMap;
    ParticleKind mapKind(Kind k,LevelIdx pLevel,KindIdx kindIdx){
        auto kVal=k.name in specialMap;
        if (kVal!is null){
            return *kVal;
        }
        char[] symb=k.symbol;
        if (symb==""){
            auto nLen=k.name.length;
            if (nLen>0){
                if (nLen>=2){
                    symb=k.name[0..2];
                    if ((symb[1]>='a' && symb[1]<='z') || (symb[1]>='A' && symb[1]<='Z')){
                        symb=symb.dup;
                    } else {
                        symb=symb[0..1];
                    }
                } else {
                    symb=symb[0..1];
                }
                if (atomFromSymbol(symb,true)is null){
                    symb=null;
                } else {
                    symb=symb.dup;
                }
            }
        }
        auto p=new ParticleKind(k.name.dup,pLevel,kindIdx,symb,k.potential.dup,((pLevel==0)?1:0),0,0);
        return p;
    }
    static ParticleKindMap defaultMap;
    static this(){
        defaultMap=new ParticleKindMap();
    }
}

ParticleSys!(T) readIn2PSys(T)(ReadSystem rIn,
    ParticleKind delegate(Kind,LevelIdx,KindIdx) pkindsMap=&ParticleKindMap.defaultMap.mapKind,
    NotificationCenter nCenter=null)
{
    sout("fullSystem & levels\n");
    // fullSystem & levels
    KindRange[] levels=new KindRange[](4);
    auto sortedPIndex=BulkArray!(PIndex)(rIn.nParticles+rIn.residui.length+rIn.chains.length+1);
    index_type[] kindStarts=new index_type[](rIn.pKinds.length+rIn.resKinds.length+rIn.chainKinds.length+2);

    size_t ii=0;
    foreach(i,pk; rIn.pKinds){
        PIndex pAtt=PIndex(cast(KindIdx)i,cast(ParticleIdx)0);
        kindStarts[cast(size_t)pAtt.kind]=ii;
        if (! pk.nextParticle.valid){
            throw new Exception("non valid nextParticle",__FILE__,__LINE__);
        }
        auto nPart=cast(size_t)pk.nextParticle.particle;
        for(size_t iPLocal=nPart;iPLocal!=0;--iPLocal){
            sortedPIndex[ii]=pAtt;
            ++ii;
            pAtt+=1;
        }
    }
    levels[0]=KindRange(cast(KindIdx)0,cast(KindIdx)rIn.pKinds.length);
    foreach(i,pk; rIn.resKinds){
        PIndex pAtt=PIndex(cast(KindIdx)(i+levels[0].kEnd),cast(ParticleIdx)0);
        kindStarts[cast(size_t)pAtt.kind]=ii;
        if (! pk.nextParticle.valid){
            throw new Exception("non valid nextParticle",__FILE__,__LINE__);
        }
        auto nPart=cast(size_t)pk.nextParticle.particle;
        for(size_t iPLocal=nPart;iPLocal!=0;--iPLocal){
            sortedPIndex[ii]=pAtt;
            ++ii;
            pAtt+=1;
        }
    }
    levels[1]=KindRange(levels[0].kEnd,cast(KindIdx)(rIn.resKinds.length+levels[0].kEnd));
    foreach(i,pk; rIn.chainKinds){
        PIndex pAtt=PIndex(cast(KindIdx)(i+levels[1].kEnd),cast(ParticleIdx)0);
        kindStarts[cast(size_t)pAtt.kind]=ii;
        if (! pk.nextParticle.valid){
            throw new Exception("non valid nextParticle",__FILE__,__LINE__);
        }
        auto nPart=cast(size_t)pk.nextParticle.particle;
        for(size_t iPLocal=nPart;iPLocal!=0;--iPLocal){
            sortedPIndex[ii]=pAtt;
            ++ii;
            pAtt+=1;
        }
    }
    levels[2]=KindRange(levels[1].kEnd,cast(KindIdx)(rIn.chainKinds.length+levels[1].kEnd));
    kindStarts[cast(size_t)levels[2].kEnd]=ii;
    sortedPIndex[ii]=PIndex(levels[2].kEnd,cast(ParticleIdx)0);// system particle
    levels[3]=KindRange(levels[2].kEnd,cast(KindIdx)(1+levels[2].kEnd));
    ++ii;
    kindStarts[cast(size_t)levels[2].kEnd+1]=ii;
    auto fullRange=KindRange(levels[0].kStart,levels[3].kEnd);
    auto fullSystem=new SubMapping("fullSystem",sortedPIndex,*cast(BulkArray!(LocalPIndex)*)cast(void*)&sortedPIndex,
        sortedPIndex,kindStarts,
        fullRange,MappingKind.Same);
    
    sout("externalOrder\n");
    // externalOrder
    auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(kindStarts[levels[0].kEnd]);
    auto lSortedPIndex=BulkArray!(PIndex)(kindStarts[levels[0].kEnd]);
    foreach(p;rIn.particles){
        gSortedLocalPIndex[kindStarts[cast(size_t)p.pIndex.kind]+cast(size_t)p.pIndex.particle]=LocalPIndex(cast(ulong)p.externalIdx);
        lSortedPIndex[p.externalIdx]=p.pIndex;
    }
    auto externalOrder=new SubMapping("externalOrder",sortedPIndex[0..kindStarts[levels[0].kEnd]],
        gSortedLocalPIndex,lSortedPIndex,kindStarts[0..1+levels[0].kEnd],KindRange(levels[0].kStart,levels[0].kEnd),
        MappingKind.Gapless);
    
    sout("particleStruct, superParticle\n");
    // particleStruct, superParticle
    index_type[] kindDims=new index_type[](cast(size_t)fullRange.kEnd);
    kindDims[]=1;
    sout("pippo1\n");
    auto particlesStruct=new SegmentedArrayStruct("particleStruct",fullSystem,fullRange,kindDims);
    auto particles=new SegmentedArray!(PIndex)(particlesStruct,sortedPIndex,fullRange,kindStarts);
    auto superParticle=new SegmentedArray!(PIndex)(particlesStruct);
    auto resShift=levels[0].kEnd;
    auto chainShift=levels[1].kEnd;
    scope nSub=new SegmentedArray!(size_t)(particlesStruct,BulkArray!(size_t).dummy,KindRange(levels[0].kEnd,levels[3].kEnd));
    nSub[]=0;
    foreach(p;rIn.particles){
        auto resIdx=PIndex(resShift+p.resIndex.kind,p.resIndex.particle);
        auto chainIdx=PIndex(chainShift+p.chainIndex.kind,p.chainIndex.particle);
        superParticle[LocalPIndex(p.pIndex),0]=resIdx;
        superParticle[LocalPIndex(resIdx),0]=chainIdx;
        *(nSub.ptrI(resIdx,0)) += 1;
    }
    foreach(p;superParticle[levels[1]].pLoop){
        atomicAdd(*(nSub.ptrI(p,0)),cast(size_t)1);
    }
    auto sysIdx=PIndex(levels[3].kStart,0);
    superParticle[levels[2]].data[]=sysIdx;
    sout("sysIdx")(sysIdx)(" levels[3].kStart")(levels[3].kStart)(" levels[3].kEnd")(levels[3].kEnd)("\n");
    sout("pippo17\n");
    *nSub.ptrI(sysIdx,0)=superParticle[levels[2]].data.length;
    sout("pippo18\n");
    
    sout("kinds particleKinds\n");
    /// kinds particleKinds
    auto kindsData=BulkArray!(ParticleKind)(rIn.pKinds.length+rIn.resKinds.length+rIn.chainKinds.length+2);
    foreach(i,pk; rIn.pKinds){
        kindsData[i]=pkindsMap(pk,cast(LevelIdx)0,pk.nextParticle.kind);
    }
    auto kindShift=levels[0].kEnd;
    foreach(i,pk; rIn.resKinds){
        kindsData[i+cast(size_t)kindShift]=pkindsMap(pk,cast(LevelIdx)1,cast(KindIdx)(kindShift+pk.nextParticle.kind));
    }
    kindShift=levels[1].kEnd;
    foreach(i,pk; rIn.chainKinds){
        kindsData[i+cast(size_t)kindShift]=pkindsMap(pk,cast(LevelIdx)2,cast(KindIdx)(kindShift+pk.nextParticle.kind));
    }
    kindsData[levels[3].kStart]=new SysKind(cast(LevelIdx)3,levels[3].kStart);
    auto kindDims2=new index_type[](kindDims.length);
    kindDims2[]=0;
    auto kindsStruct=new SegmentedArrayStruct("kindsStruct",fullSystem,fullRange,kindDims2,SegmentedArrayStruct.Flags.Min1);
    auto particleKinds=new SegmentedArray!(ParticleKind)(kindsStruct,kindsData);
    
    sout("subparticles\n");
    /// subparticles
    index_type[] nSubparticles=new index_type[](cast(size_t)levels[3].kEnd);
    nSubparticles[]=index_type.max;
    foreach(lIdx,val;nSub.pLoop){
        auto kind=cast(size_t)lIdx.kind;
        auto nP=nSubparticles[kind];
        if (nP==index_type.max){
            nP=atomicCAS(nSubparticles[kind],cast(index_type)val,index_type.max);
        }
        if (nP!=val && nP!=index_type.max){
            throw new Exception(collectAppender(delegate void(CharSink sink){
                auto s=dumperP(sink);
                s("Found different number of subparticles:")(nP)(" vs ")(val)
                (" for particle of kind ")(kindsData[kind])("\n");
            }),__FILE__,__LINE__);
        }
    }
    foreach (i,v;nSubparticles){
        if (v!=index_type.max) {
            auto oldV=kindsData[i].subParticles;
            if(oldV!=0 && oldV!=size_t.max && oldV!=v && oldV!=index_type.max){
                throw new Exception(collectAppender(delegate void(CharSink sink){
                    dumperP(sink)("subParticles was set to ")(oldV)(" but found ")(v)
                        (" subparticles for particles of kind ")(kindsData[i])("\n");
                }),__FILE__,__LINE__);
            }
            kindsData[i].subParticles=v;
        }
    }

    auto subParticlesStruct=new SegmentedArrayStruct("subparticleStruct",fullSystem,KindRange(levels[1].kStart,levels[3].kEnd),nSubparticles);
    auto subParticleIdxs=new SegmentedArray!(PIndex)(subParticlesStruct);
    nSub[]=0;
    foreach(lIdx,superP;superParticle.sLoop){ // sequential, we need to guarantee a deterministic result
        auto idxAtt=nSub.ptrI(lIdx,0);
        *(subParticleIdxs.ptrI(superP,*idxAtt))=PIndex(lIdx);
        ++(*idxAtt);
    }
    Exception e;
    foreach(lIdx,nPart;nSub.pLoop){
        if (nSubparticles[cast(size_t)lIdx.kind]!=cast(index_type)nPart){
            e=new Exception(collectAppender(delegate void(CharSink sink){
                dumperP(sink)("internal error inconsistent number of subparticles:")
                    (nSubparticles[cast(size_t)lIdx.kind])("vs")(nPart)("\n");
            }),__FILE__,__LINE__);
            break;
        }
    }
    if (e!is null) throw e;
    if (nSub.data.guard!is null) nSub.data.guard.release();
    
    sout("sysStruct\n");
    // sysStruct
    auto sysStruct=new SysStruct(rIn.name, fullSystem, externalOrder, levels,
         particlesStruct, particles, superParticle, subParticlesStruct,
         subParticleIdxs, kindsStruct, particleKinds);
    
    DynamicsVars!(T) dynVars;
    
    if (nCenter is null) nCenter=new NotificationCenter();
    
    sout("particleSystem\n");
    // particleSystem
    auto pSys=new ParticleSys!(T)(0,rIn.name,sysStruct,nCenter);
    
    pSys.reallocStructs();
    pSys.sysStructChanged();
    
    pSys.checkX();
    Matrix!(T,3,3) h;
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            h[i,j]=scalar!(T)(rIn.cell[i][j]);
        }
    }
    pSys.dynVars.cell=new Cell!(T)(h,rIn.periodic,Vector!(T,3)(rIn.x0[0],rIn.x0[1],rIn.x0[2]));
    
    sout("setPos\n");
    auto posV=pSys.dynVars.pos;
    foreach (p;rIn.particles){
        posV[LocalPIndex(p.pIndex),0]=Vector!(T,3)(p.pos[0],p.pos[1],p.pos[2]);
    }
    
    sout("cellUpdateNotification\n");
    pSys.cellChanged();
    sout("posUpdateNotification\n");
    pSys.positionsChanged();
    sout("pSys setup end\n");
    
    return pSys;
}