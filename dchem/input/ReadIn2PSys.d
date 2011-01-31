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

class ParticleKindMap{
    ParticleKind[char[]] specialMap;
    ParticleKind mapKind(Kind k,LevelIdx pLevel,KindIdx kindIdx){
        auto kVal=k.name in specialMap;
        if (kVal!is null){
            return *kVal;
        }
        char[] symb=k.symbol;
        if (symb.length==0){
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
                    symb=k.name[0..1];
                }
                if (atomFromSymbol(symb,true)is null){
                    symb=null;
                } else {
                    symb=symb.dup;
                }
            }
        }
        auto p=new ParticleKind(k.name.dup,pLevel,kindIdx,symb,k.potential.dup);
        index_type elIdx,delIdx;
        if (pLevel==0) p.posEls.addElements(1,1,false,elIdx,delIdx); // add a position
        return p;
    }
    static ParticleKindMap defaultMap;
    static this(){
        defaultMap=new ParticleKindMap();
    }
}

ParticleSys!(T) readIn2PSys(T)(ReadSystem rIn,
    ParticleKind delegate(Kind,LevelIdx,KindIdx) pkindsMap=null,
    NotificationCenter nCenter=null)
{
    if (pkindsMap is null) pkindsMap=&ParticleKindMap.defaultMap.mapKind;
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
    
    // externalOrder
    auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(kindStarts[levels[0].kEnd]);
    auto lSortedPIndex=BulkArray!(PIndex)(kindStarts[levels[0].kEnd]);
    foreach(p;rIn.particles){
        lSortedPIndex[kindStarts[cast(size_t)p.pIndex.kind]+cast(size_t)p.pIndex.particle]=PIndex(cast(ulong)p.externalIdx);
        gSortedLocalPIndex[p.externalIdx]=LocalPIndex(p.pIndex);
    }
    auto sortedPIndex2=BulkArray!(PIndex)(kindStarts[levels[0].kEnd]);
    foreach (i,ref v;sortedPIndex2){
        v=PIndex(0,i);
    }
    auto externalOrder=new SubMapping("externalOrder",sortedPIndex2,
        gSortedLocalPIndex,lSortedPIndex,kindStarts[0..1+levels[0].kEnd],levels[0],
        MappingKind.Gapless);
    
    // particleStruct, superParticle
    index_type[] kindDims=new index_type[](cast(size_t)fullRange.kEnd);
    kindDims[]=1;
    auto particlesStruct=new SegmentedArrayStruct("particlesStruct",fullSystem,fullRange,kindDims);
    auto pMapIdx=new SegArrMemMap!(PIndex)(particlesStruct,KindRange.all);
    auto particles=pMapIdx.newArray(sortedPIndex,fullRange);
    auto superParticle=pMapIdx.newArray();
    auto resShift=levels[0].kEnd;
    auto chainShift=levels[1].kEnd;
    auto pMapKSizeT=new SegArrMemMap!(size_t)(particlesStruct,KindRange(levels[0].kEnd,levels[3].kEnd));
    scope nSub=pMapKSizeT.newArray();
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
    superParticle[levels[2]].support[]=sysIdx;
    *nSub.ptrI(sysIdx,0)=superParticle[levels[2]].support.length;
    
    /// kinds particleKinds
    auto kindsData=BulkArray!(ParticleKind)(rIn.pKinds.length+rIn.resKinds.length+rIn.chainKinds.length+1);
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
    auto kindMapPKind=new SegArrMemMap!(ParticleKind)(kindsStruct,KindRange.all);
    auto particleKinds=kindMapPKind.newArray(kindsData);
    
    /// subparticles
    index_type[] nSubparticles=new index_type[](cast(size_t)levels[3].kEnd);
    nSubparticles[]=index_type.max;
    foreach(inP,pIdx,lIdx,val;nSub.pLoop){
        auto kind=cast(size_t)lIdx.kind;
        auto nP=nSubparticles[kind];
        if (nP==index_type.max){
            nP=atomicCAS(nSubparticles[kind],cast(index_type)val,index_type.max);
        }
        if (nP!=val && nP!=index_type.max){
            throw new Exception(collectAppender(delegate void(CharSink sink){
                auto s=dumper(sink);
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
                    dumper(sink)("subParticles was set to ")(oldV)(" but found ")(v)
                        (" subparticles for particles of kind ")(kindsData[i])("\n");
                }),__FILE__,__LINE__);
            }
            kindsData[i].subParticles=v;
        }
    }

    auto subParticlesStruct=new SegmentedArrayStruct("subParticlesStruct",fullSystem,KindRange(levels[1].kStart,levels[3].kEnd),nSubparticles[levels[1].kStart..levels[3].kEnd]);
    auto subPMapPIndex=new SegArrMemMap!(PIndex)(subParticlesStruct);
    auto subParticleIdxs=subPMapPIndex.newArray();
    nSub[]=0;
    foreach(inP,pIdx,lIdx,superP;superParticle[KindRange(levels[0].kStart,levels[2].kEnd)].sLoop){ // sequential, we need to guarantee a deterministic result
        nSub.dtype* idxAtt;
        idxAtt=nSub.ptrI(superP,0);
        *(subParticleIdxs.ptrI(superP,*idxAtt))=PIndex(lIdx);
        ++(*idxAtt);
    }
    Exception e;
    foreach(inP,pIdx,lIdx,nPart;nSub.pLoop){
        if (nSubparticles[cast(size_t)lIdx.kind]!=cast(index_type)nPart){
            e=new Exception(collectAppender(delegate void(CharSink sink){
                dumper(sink)("internal error inconsistent number of subparticles:")
                    (nSubparticles[cast(size_t)lIdx.kind])("vs")(nPart)("\n");
            }),__FILE__,__LINE__);
            break;
        }
    }
    if (e!is null) throw e;
    pMapKSizeT.rmUser();
    nSub.giveBack();
    
    // sysStruct
    auto sysStruct=new SysStruct(rIn.name, fullSystem, externalOrder, levels,
         particlesStruct, particles, superParticle, subParticlesStruct,
         subParticleIdxs, kindsStruct, particleKinds);
    
    if (nCenter is null) nCenter=new NotificationCenter();
    
    // particleSystem
    auto pSys=new ParticleSys!(T)(0,rIn.name,sysStruct,nCenter);
    pSys.pKindsInitialSetup(); // first setup of particle kinds
    pSys.sysStructChanged();
    
    pSys.checkX();
    Matrix!(T,3,3) h;
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            h[i,j]=scalar!(T)(rIn.cell[i][j]);
        }
    }
    pSys.dynVars.x.cell=new Cell!(T)(h,rIn.periodic,Vector!(T,3)(rIn.x0[0],rIn.x0[1],rIn.x0[2]));
    
    auto posV=pSys.dynVars.x.pos;
    foreach (p;rIn.particles){
        posV[LocalPIndex(p.pIndex),0]=Vector!(T,3)(p.pos[0],p.pos[1],p.pos[2]);
    }
    
    pSys.cellChanged();
    pSys.positionsChanged();
    
    return pSys;
}


/// builds an artificial particle system without particles, only with systemwide attributes
ParticleSys!(T) artificialPSys(T)(size_t nPos, size_t nOrient,size_t nDof,
    NotificationCenter nCenter=null)
{
    // fullSystem & levels
    KindRange[] levels=new KindRange[](4);
    auto sortedPIndex=BulkArray!(PIndex)(1);
    index_type[] kindStarts=new index_type[](2);

    kindStarts=[cast(index_type)0,1];
    for (size_t ipart=0;ipart<3;++ipart){
        levels[ipart]=KindRange(cast(KindIdx)0,cast(KindIdx)0);
    }
    levels[3]=KindRange(cast(KindIdx)0,cast(KindIdx)1);
    sortedPIndex[0]=PIndex(3,0); // the system particle
    
    auto fullRange=KindRange(levels[0].kStart,levels[$-1].kEnd);
    auto fullSystem=new SubMapping("fullSystem",sortedPIndex,*cast(BulkArray!(LocalPIndex)*)cast(void*)&sortedPIndex,
        sortedPIndex,kindStarts,
        fullRange,MappingKind.Same);
    
    // externalOrder
    auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(0);
    auto lSortedPIndex=BulkArray!(PIndex)(0);

    auto externalOrder=new SubMapping("externalOrder",sortedPIndex[0..kindStarts[levels[0].kEnd]],
        gSortedLocalPIndex,lSortedPIndex,kindStarts[0..1+levels[0].kEnd],levels[0],
        MappingKind.Gapless);
    
    // particleStruct, superParticle
    index_type[] kindDims=new index_type[](cast(size_t)fullRange.kEnd);
    kindDims[]=1;
    auto particlesStruct=new SegmentedArrayStruct("particlesStruct",fullSystem,fullRange,kindDims);
    auto partMapPIndex=new SegArrMemMap!(PIndex)(particlesStruct);
    auto particles=partMapPIndex.newArray(sortedPIndex);
    auto superParticle=partMapPIndex.newArray();
    
    superParticle[LocalPIndex(0,0),0]=PIndex();

    /// kinds particleKinds
    auto kindsData=BulkArray!(ParticleKind)(1);
    auto p=new SysKind(cast(LevelIdx)3,levels[3].kStart);
    
    kindsData[0]=p;
    index_type elIdx,delIdx;
    p.posEls.addElements(nPos,nPos,false,elIdx,delIdx);
    p.orientEls.addElements(nOrient,0,true,elIdx,delIdx);
    p.dofEls.addElements(nDof,nDof,false,elIdx,delIdx);
    p.dofEls.addElements(0,3*nOrient,true,elIdx,delIdx);
    auto kindDims2=new index_type[](kindDims.length);
    kindDims2[]=0;
    auto kindsStruct=new SegmentedArrayStruct("kindsStruct",fullSystem,fullRange,kindDims2,SegmentedArrayStruct.Flags.Min1);
    auto kindsMapPKind=new SegArrMemMap!(ParticleKind)(kindsStruct);
    auto particleKinds=kindsMapPKind.newArray(kindsData);
    
    /// subparticles
    index_type[] nSubparticles=new index_type[](1);
    nSubparticles[0]=0;
    foreach (i,v;nSubparticles){
        if (v!=index_type.max) {
            auto oldV=kindsData[i].subParticles;
            if(oldV!=0 && oldV!=size_t.max && oldV!=v && oldV!=index_type.max){
                throw new Exception(collectAppender(delegate void(CharSink sink){
                    dumper(sink)("subParticles was set to ")(oldV)(" but found ")(v)
                        (" subparticles for particles of kind ")(kindsData[i])("\n");
                }),__FILE__,__LINE__);
            }
            kindsData[i].subParticles=v;
        }
    }

    auto subParticlesStruct=new SegmentedArrayStruct("subParticlesStruct",fullSystem,KindRange(levels[1].kStart,levels[$-1].kEnd),nSubparticles[levels[1].kStart..levels[$-1].kEnd]);
    auto subPMapPIdx=new SegArrMemMap!(PIndex)(subParticlesStruct);
    auto subParticleIdxs=subPMapPIdx.newArray();
    auto pMapSizeT=new SegArrMemMap!(size_t)(particlesStruct,KindRange(levels[0].kEnd,levels[$-1].kEnd));
    scope nSub=pMapSizeT.newArray();
    nSub[]=0;
    foreach(inP,pIdx,lIdx,superP;superParticle[KindRange(levels[0].kStart,levels[$-2].kEnd)].sLoop){ // sequential, we need to guarantee a deterministic result
        nSub.dtype* idxAtt;
        idxAtt=nSub.ptrI(superP,0);
        *(subParticleIdxs.ptrI(superP,*idxAtt))=PIndex(lIdx);
        ++(*idxAtt);
    }
    Exception e;
    foreach(inP,pIdx,lIdx,nPart;nSub.pLoop){
        if (nSubparticles[cast(size_t)lIdx.kind]!=cast(index_type)nPart){
            e=new Exception(collectAppender(delegate void(CharSink sink){
                dumper(sink)("internal error inconsistent number of subparticles:")
                    (nSubparticles[cast(size_t)lIdx.kind])("vs")(nPart)("\n");
            }),__FILE__,__LINE__);
            break;
        }
    }
    if (e!is null) throw e;
    pMapSizeT.rmUser();
    nSub.giveBack(); nSub=null;
    
    // sysStruct
    auto sysStruct=new SysStruct("sys", fullSystem, externalOrder, levels,
         particlesStruct, particles, superParticle, subParticlesStruct,
         subParticleIdxs, kindsStruct, particleKinds);
    
    if (nCenter is null) nCenter=new NotificationCenter();
    
    // particleSystem
    auto pSys=new ParticleSys!(T)(0,"sys",sysStruct,nCenter);
    pSys.pKindsInitialSetup(); // first setup of particle kinds
    pSys.sysStructChanged();
    
    pSys.checkX();
    Matrix!(T,3,3) h=Matrix!(T,3,3).identity;

    pSys.dynVars.x.cell=new Cell!(T)(h,[0,0,0],Vector!(T,3).zero);
    
    auto posV=pSys.dynVars.x;
    posV.nullCell[]=0;
    
    pSys.cellChanged();
    pSys.positionsChanged();
    
    return pSys;
}