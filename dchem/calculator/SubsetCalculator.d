/// builds a calculator that has only a subset as active part
module dchem.calculator.SubsetCalculator;
import dchem.sys.ParticleRange;
import dchem.Common;
import blip.serialization.Serialization;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import blip.BasicModels;
import blip.container.BulkArray;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import dchem.input.RootInput;
import dchem.sys.DynVars;

/// selection mapping of a particle
struct ParticleDynvarsSelect{
    enum Selection:uint{
        x=1,
        y=1<<1,
        z=1<<2,
        orient=1<<3,
        dof=1<<4,
        dofIdxs=1<<5,
        xyz=x|y|z,
        direct=x|y|z|orient|dof,
    }
    uint selection;
    ulong[] dofIdxs;
    size_t posDofs,posDDos;
    mixin(serializeSome("dchem.ParticleDynvarsSelect",`selection|dofIdxs|posDofs`));
    mixin printOut!();
}

/// represent a given subset of the dynamic variables of a system
struct SubsetSelection{
    bool x=true,y=true,z=true;
    bool orient=true;
    bool dof=true;
    ulong[] dofRanges;
    ulong[] dofIdxs;
    
    InputField[] particleRanges;
    
    mixin(serializeSome("dchem.SubsetSelection",`
    x: select x component (true)
    y: select y component (true)
    z: select z component (true)
    orient: select the orientation (true)
    dof: select the degrees of freedom (true)
    dofIdxs: indexes of degrees of freedom to select (if empty the full range is assumed)
    particleRanges: the ranges of particles that have this selection
    `));
    mixin printOut!();
    
    bool verify(CharSink s,char[] fieldName){
        bool res=true;
        foreach (p;particleRanges){
            if (p is null){
                dumper(s)("particleRanges should not be null in field ")(fieldName);
                res=false;
            } else if ((cast(ParticleRange)p)is null){
                dumper(s)("particleRanges should be of type ParticleRange in field ")(fieldName);
                res=false;
            }
        }
        return res;
    }
    
    ParticleKind newParticleKind(ParticleKind oldPKind,KindIdx newKIdx,LevelIdx newLevel,char[]newName,out ParticleDynvarsSelect pSelect){
        auto p=new ParticleKind(newName,newLevel,newKIdx,oldPKind.symbol,oldPKind.potential);
        alias ParticleDynvarsSelect.Selection Sel;
        if (x) pSelect.selection|=Sel.x;
        if (y) pSelect.selection|=Sel.y;
        if (z) pSelect.selection|=Sel.z;
        if (orient){
            pSelect.selection|=Sel.orient;
            p.orientEls=oldPKind.orientEls;
        }
        if (dof){
            pSelect.selection|=Sel.dof;
            auto ndof=oldPKind.dofEls.nElements;
            auto nddof=oldPKind.dofEls.nDelements;
            if (dofIdxs.length!=0){
                {
                    scope visited=new bool[](ndof);
                    visited[]=false;
                    foreach(d;dofIdxs){
                        assert(d>=0&&d<np);
                        assert(!visited[d]);
                        visited[d]=true;
                    }
                }
                pSelect.selection|=Sel.dofIdxs
                if (oldPKind.posEls.derivMap!=ParticleKind.DerivMap.SimpleMap){
                    new Exception("SubsetSelection doesn't know how to handle complex map of dofs with selection in particle "
                        ~oldPKind.name,__FILE__,__LINE__);
                }
                pSelect.dofIdxs=dofIdxs;
            }
        }
        {
            auto np=oldPKind.posEls.nElements;
            auto dnp=oldPKind.posEls.nDelements;
            if ((x&&y&&z)||(pSelect.selections&Sel.xyz)==0)){
                if (oldPKind.posEls.derivMap!=ParticleKind.DerivMap.SimpleMap){
                    new Exception("SubsetSelection doesn't know how to handle complex map in pos in particle "
                        ~oldPKind.name,__FILE__,__LINE__);
                }
                if (np!=dnp){
                    new Exception("SubsetSelection doesn't know how to handle nElements!=nDelements in pos in particle "
                        ~oldPKind.name,__FILE__,__LINE__);
                }
            }
            size_t elIdx,delIdx
            switch (pSelect.selections&Sel.xyz){
                case Sel.xyz:
                    p.posEls=oldPKind.posEls;
                    break;
                case (Sel.x|Sel.y),(Sel.x|Sel.z),(Sel.y|Sel.z):
                    p.dofEls.addElements(2*np,2*np,false,pSelect.posDofs,pSelect.posDDofs);
                    break;
                case Sel.x,Sel.y,Sel.z:
                    p.dofEls.addElements(2*np,2*np,false,pSelect.posDofs,pSelect.posDDofs);
                    break;
                case 0:
                    break
                default:
                    assert(0);
            }
        }
        
    }
}

/// returns a calculator that uses a subset of another calculator
class SubsetCalculator: Method{
    SubsetSelection[] selections;
    InputField calculator;
    bool trimUnspecifiedSuperParticles=true;
    bool singleSys=true;
    
    SysStruct fullSysStruct; /// structure of the full system (if singleSys is true)
    SysStruct sysStruct; /// structure of this subsetted system (if singleSys is true)
    SubMapping subSysMap; /// map particles in the whole system (PIndex) to/from this (LocalPIndex) (if singleSys is true)
    BulkArray!(ParticleDynvarsSelect) dynVarMap; /// mapping of the dynVarMap for each (local) kind (if singleSys is true)
    RLock gLock;
    
    this(){
        gLock=new RLock();
    }
    
    Method subMethod(){
        auto m=cast(Method)calculator;
        assert(m!is null);
        return m;
    }
    mixin(serializeSome("dchem.SubsetCalculator",`
    selections: the selected degrees of freedom
    calculator: the calculator to subset
    trimUnspecifiedSuperParticles: if unspecified superparticle should have their degrees of freedom removed (true)
    singleSys: if there is a single system for all the contexts, and thus the mapping and new sysStruct can be cached (true)`));
    mixin printOut!();
    mixin myFieldMixin!();
    /// activates the method in the context of the given parallel environment
    /// this can be used to do setups shared by all contexts of this method
    /// it might be called several times, but should always have the same arguments
    void setup(LinearComm pEnv, CharSink log){
        subMethod.setup(pEnv,log);
    }
    /// gets a calculator to perform calculations with this method, if possible reusing the given history
    /// if wait is true waits until a context is available
    CalculationContext getCalculator(bool wait,ubyte[]history){
        auto cCtx=subMethod.getCalculator(wait,history);
        auto prec=cCtx.activePrecision;
        switch(prec){
            case Precision.Real:
                return new SubsetContext!(Real)(this,cCtx);
            case Precision.LowP:
                return new SubsetContext!(LowP)(this,cCtx);
            default:
                assert(0,"unexpected precision");
        }
    }
    /// drops the history with the given id
    void dropHistory(ubyte[]history){
        subMethod.dropHistory(history);
    }
    /// clears all history
    void clearHistory(){
        subMethod.clearHistory();
    }
    
    bool verify(CharSink s){
        bool res=true;
        foreach (s;selections){
            res=res&&s.verify(s,myFieldName);
        }
        if (calculator is null){
            dumper(s)("calculator should not be null in field ")(myFieldName);
            res=false;
        } else if ((cast(Calculator)p)is null){
            dumper(s)("calculator should be a Calculator in field ")(myFieldName);
            res=false;
        }
        return res;
    }
}

class SubsetContext(T):CalcContext{
    SubsetCalculator input;
    CalculationContext subContext;
    SysStruct sysStruct; /// structure of this subsetted system
    SubMapping subSysMap; /// map particles in the whole system (PIndex) to/from this (LocalPIndex)
    BulkArray!(ParticleDynvarsSelect) dynVarMap; /// mapping of the dynVarMap for each (local) kind
    
    static struct PK{
        size_t isub=size_t.max;
        KindIdx pKind;
        KindIdx pKindNew;
        static PK opCall(size_t isub,KindIdx pKind,KindIdx pKindNew=KindIdx.init){
            PK res;
            res.isub=isub;
            res.pKind=pKind;
            res.pKindNew=pKindNew;
            return res;
        }
        mixin(serializeSome("schem.SubsetContext.PK","isub|pKind|id|firstPK"));
        mixin printOut!();
    }
    static struct PKV{
        PK pk;
        PIndex[] subP;
        size_t nParticles;
        mixin(serializeSome("schem.SubsetContext.PK","pk|subP"));
        mixin printOut!();
    }
    
    this(SubsetCalculator input,CalculationContext subContext,NotificationCenter nCenter=null){
        this.input=input;
        this.subContext=subContext;
        
        auto fullSysStruct=subContext.sysStruct;
        
        if (input.singleSys){
            input.gLock.lock();
            scope(exit){
                input.gLock.unlock();
            }
            if (input.sysStruct is null){
                calcSysStruct(fullSysStruct);
                input.sysStruct=sysStruct;
                input.subSysMap=subSysMap;
                input.dynVarMap=dynVarMap;
                input.fullSysStruct=fullSysStruct;
            } else {
                sysStruct=input.sysStruct;
                subSysMap=input.subSysMap;
                dynVarMap=input.dynVarMap;
                // skip check?
                assert(fullSysStruct==input.fullSysStruct,"full struct mismatch");
            }
        } else {
            calcSysStruct();
        }
        
        if (nCenter is null) nCenter=new NotificationCenter();

        // particleSystem
        auto pSys=new ParticleSys!(T)(0,,sysStruct,nCenter);
        pSys.pKindsInitialSetup(); // first setup of particle kinds
        pSys.sysStructChanged();

        pSys.checkX();
        Matrix!(T,3,3) h;
        for (int i=0;i<3;++i){
            for (int j=0;j<3;++j){
                h[i,j]=scalar!(T)(rIn.cell[i][j]);
            }
        }
        pSys.dynVars.x.cell=subContext.ref;

        auto posV=pSys.dynVars.x.pos;
        foreach (p;rIn.particles){
            posV[LocalPIndex(p.pIndex),0]=Vector!(T,3)(p.pos[0],p.pos[1],p.pos[2]);
        }

        pSys.cellChanged();
        pSys.positionsChanged();

    }
    
    /// calculates sysStruct, subSysMap, and dynVarMap using the given subSys as full system
    void calcSysStruct(SysStruct subSys){
        auto pKinds=new SegmentedArray!(PK)(subSys.particlesStruct);
        pKinds[]=PK.init;
        foreach(isub,subset;input.selections){
            foreach(pIdx;subset){
                if (pKinds[pIdx]!=PK.init){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("double use of ")(pIdx)(" from subset ")(pKinds[pIdx].isub)
                            (" and ")(isub);
                    }),__FILE__,__LINE__);
                }
                pKinds[pIdx]=PK(isub,pIdx.kind);
            }
        }
        auto sLevels=subSys.levels;
        KindIdx nextK=0;
        auto newKinds=new HashMap!(PK,PKV[])();
        size_t nParticles;
        foreach(pIdx,lIdx,ref pk;pKinds[sLevels[0]]){
            if (pk.pKind!=KindIdx.init){
                ++nParticles;
                auto nK=pk in newKinds;
                if (nK is null){
                    auto res=new PKV[](1);
                    res[0].pk=pk;
                    res[0].pk.pKindNew=nextK++;
                    ++(res[0].nParticles);
                    newKinds[pk]=res;
                    pk.pKindNew=res[0].pk.pKindNew;
                } else {
                    ++((*nK)[0].nParticles);
                    pk.pKindNew=(*nK)[0].pk.pKindNew;
                }
            }
        }
        KindRange[] levels=new KindRange[](sLevels.length);
        levels[0].kStart=0;
        levels[0].kEnd=nextK;
        for (int ilevel=1;ilevel<sLevels.length;++ilevel){
            levels[ilevel].kStart=nextK;
            foreach(pIdx,lIdx,ref pk;pKinds[sLevels[0]]){
                auto subNow=subSys.subParticles[pIdx];
                bool hasSome=false;
                foreach (pIdx;subNow){
                    if (pKinds[pIdx].pKind!=KindIdx.init){
                        hasSome=true;
                        break;
                    }
                }
                if (!hasSome) continue;
                ++nParticles;
                if (pk.pKind==KindIdx.init) pk.pKind=pIdx;
                auto nK=pk in newKinds;
                if (nK !is null){
                    bool same=false;
                    foreach(ref pkv;*nK){
                        same=true;
                        foreach(iSub,subPk,pkv.subP){
                            if (pKinds[subPk]!=pKinds[subNow[iSub]]){
                                same=false;
                                break;
                            }
                        }
                        if (same){
                            ++pkv.nParticles;
                            pk.pKindNew=pkv.pk.pKindNew;
                            break;
                        }
                    }
                    if (!same){
                        PKV res;
                        res.pk=pk;
                        res.pk.pKindNew=nextK++;
                        res.subP=subNow;
                        ++res.nParticles;
                        (*nK)~=res;
                        pk.pKindNew=res.pk.pKindNew;
                    }
                } else {
                    auto res=new PKV[](1);
                    res[0].pk=pk;
                    res[0].pk.pKindNew=nextK++;
                    res[0].subP=subNow;
                    ++res[0].nParticles;
                    newKinds[pk]=res;
                    pk.pKindNew=res[0].pk.pKindNew;
                }
            }
            levels[ilevel].kEnd=nextK;
        }
        /// lookup table for pkvs
        scope pkvLookup= new PKV[](cast(size_t)nextK);
        foreach(pk,pkvs;newKinds){ // worth building a lokup table?
            foreach(pkv;pkvs){
                pkvLookup[cast(size_t)pkv.pKindNew]=pkv;
            }
        }
        
        // set up sysStruct stuff
        
        auto sortedPIndex=BulkArray!(PIndex)(nParticles);
        index_type[] kindStarts=new index_type[](cast(size_t)nextK+1);
        size_t ii=0;
        for(size_t i;i<cast(size_t)nextK;++i){
            PIndex pAtt=PIndex(cast(KindIdx)i,cast(ParticleIdx)0);
            kindStarts[i]=ii;
            auto pkvAtt=pkvLookup[i];
            for(size_t iPLocal=pkvAtt.nParticles;iPLocal!=0;--iPLocal){
                sortedPIndex[ii]=pAtt;
                ++ii;
                pAtt+=1;
            }
        }
        kindStarts[$-1]=ii;
        assert(ii==nParticles);

        auto fullRange=KindRange(levels[0].kStart,levels[$-1].kEnd);
        auto fullSystem=new SubMapping("fullSystem",sortedPIndex,*cast(BulkArray!(LocalPIndex)*)cast(void*)&sortedPIndex,
            sortedPIndex,kindStarts,
            fullRange,MappingKind.Same);

        // mapping to subSys (subSys is PIndex, this PIndex is LocalPIndex)
        {
            size_t iPIdx=0;
            BulkArray!(PIndex) sortedPIndex2=BulkArray!(PIndex)(nParticles);
            auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(nParticles);
            auto lSortedPIndex=BulkArray!(PIndex)(nParticles);
            auto nParts=new size_t[](cast(size_t)nextK);
            foreach(pIdx,lIdx,pk;pKinds[sLevels[0]]){
                auto newK=pk.pKindNew;
                if (newK!=KindIdx.init){
                    sortedPIndex2[iPIdx]=pIdx;
                    size_t idx=nParts[cast(size_t)newK]++;
                    size_t lPos=kindStarts[cast(size_t)newK]+idx;
                    gSortedLocalPIndex[iPIdx]=LocalPIndex(newK,cast(ParticleIdx)idx);
                    lSortedPIndex[lPos]=pIdx;
                    ++iPIdx;
                }
            }
            assert(iPIdx==nParticles);
            subSysMap=new SubMapping("subsysMap",sortedPIndex2,
                gSortedLocalPIndex,lSortedPIndex,kindStarts,fullRange,MappingKind.Generic);
        }
        
        // externalOrder, refers to the real external order, use subSysMap???
        SubMapping externalOrder;
        {
            auto sExternalOrder=subSys.externalOrder;
            auto sortedPIndex2=BulkArray!(PIndex)(kindStarts[levels[0].kEnd]);
            auto gSortedLocalPIndex=BulkArray!(LocalPIndex)(kindStarts[levels[0].kEnd]);
            auto lSortedPIndex=BulkArray!(PIndex)(kindStarts[levels[0].kEnd]);
            auto nParts=new size_t[](cast(size_t)nextK);
            size_t iPIdx=0;
            foreach(pIdx,lIdx,pk;pKinds[sLevels[0]]){
                auto newK=pk.pKindNew;
                if (newK!=KindIdx.init){
                    size_t idx=nParts[cast(size_t)newK]++;
                    PIndex extIdx=sExternalOrder[LocalPIndex(pIdx)];
                    size_t lPos=kindStarts[cast(size_t)newK]+idx;
                    sortedPIndex2[iPIdx]=extIdx;
                    lSortedPIndex[lPos]=extIdx;
                    gSortedLocalPIndex[iPIdx]=LocalPIndex(newK,cast(ParticleIdx)idx);
                }
            }
            externalOrder=new SubMapping("externalOrder",sortedPIndex[0..kindStarts[levels[0].kEnd]],
                gSortedLocalPIndex,lSortedPIndex,kindStarts[0..1+levels[0].kEnd],KindRange(levels[0].kStart,levels[0].kEnd),
                MappingKind.Generic);
        }

        // particleStruct, superParticle
        index_type[] kindDims=new index_type[](cast(size_t)fullRange.kEnd);
        kindDims[]=1;
        auto particlesStruct=new SegmentedArrayStruct("particlesStruct",fullSystem,fullRange,kindDims);
        auto pMapIdx=new SegArrMemMap!(PIndex)(particlesStruct,KindRange.all);
        auto particles=pMapIdx.newArray(sortedPIndex,fullRange);
        auto superParticle=pMapIdx.newArray();
        auto pMapKSizeT=new SegArrMemMap!(size_t)(particlesStruct,KindRange(levels[0].kEnd,levels[$-1].kEnd));
        scope nSub=pMapKSizeT.newArray();
        nSub[]=0;
        foreach(pIdx,lIdx,ref superP;pMapIdx){
            superP=PIndex(subSysMap[subSys.superParticle[subSysMap[lIdx]]]); // probably it would be better to loop on subSys.superParticle...
            *(nSub.ptrI(superP,0)) += 1;
        }

        /// kinds particleKinds, dynVarMap
        auto kindsData=BulkArray!(ParticleKind)(cast(size_t)nextK);
        dynVarMap=BulkArray!(ParticleDynvarsSelect)(cast(size_t)nextK);
        SubsetSelection defaultSelection;
        if (trimUnspecifiedSuperParticles){
            defaultSelection.x=false;
            defaultSelection.y=false;
            defaultSelection.z=false;
            defaultSelection.orient=false;
            defaultSelection.dof=false;
        }
        auto names=new Set!(string)();
        foreach(pkv; pkvLookup){
            SubsetSelection *selection;
            if (pkv.isub==size_t.max){
                selection=&defaultSelection;
            } else {
                selection=&(selections[pkv.isub]);
            }
            ParticleKind oldPKind=subSys.kinds[cast(size_t)pkv.pKind];
            string newName=oldPKind.name;
            {
                string baseName=newName;
                char[128] buf;
                auto arr=lGrowableArray(buf,0,GASharing.Local);
                for (size_t iName=0;true;++iName){
                    auto n=newName in names;
                    if (n is null){
                        newName=newName.dup
                        names.add(newName);
                        break;
                    }
                    arr.clearData();
                    dumper(&arr.appendArr)(baseName)("_")(iName);
                    newName=arr.data;
                }
            }
            auto p=selection.newParticleKind(oldPKind,pkv.pKindNew,oldPKind.level,newName,
                *dynVarMap.ptrI(cast(size_t)pkv.pKindNew);
            kindsData[cast(size_t)pkv.pKindNew]=p;
        }
        auto kindDims2=new index_type[](kindDims.length);
        kindDims2[]=0;
        auto kindsStruct=new SegmentedArrayStruct("kindsStruct",fullSystem,fullRange,kindDims2,SegmentedArrayStruct.Flags.Min1);
        auto kindMapPKind=new SegArrMemMap!(ParticleKind)(kindsStruct,KindRange.all);
        auto particleKinds=kindMapPKind.newArray(kindsData);

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

        auto subParticlesStruct=new SegmentedArrayStruct("subParticlesStruct",fullSystem,KindRange(levels[1].kStart,levels[$-1].kEnd),nSubparticles[levels[1].kStart..levels[$-1].kEnd]);
        auto subPMapPIndex=new SegArrMemMap!(PIndex)(subParticlesStruct);
        auto subParticleIdxs=subPMapPIndex.newArray();
        nSub[]=0;
        foreach(lIdx,superP;superParticle[KindRange(levels[0].kStart,levels[2].kEnd)].sLoop){ // sequential, we need to guarantee a deterministic result
            nSub.dtype* idxAtt;
            idxAtt=nSub.ptrI(superP,0);
            *(subParticleIdxs.ptrI(superP,*idxAtt))=PIndex(lIdx);
            ++(*idxAtt);
        }
        Exception e;
        foreach(lIdx,nPart;nSub.pLoop){
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
        sysStruct=new SysStruct(rIn.name, fullSystem, externalOrder, levels,
             particlesStruct, particles, superParticle, subParticlesStruct,
             subParticleIdxs, kindsStruct, particleKinds);
    }

    ConstraintI!(Real) constraintsReal(){
        auto subC=subContext.constraintsReal();
        if (subC !is null){
            throw new Exception("constraints in the full system not yet supported",__FILE__,__LINE__);
        }
        return null;
    }
    ConstraintI!(LowP) constraintsLowP(){
        auto subC=subContext.constraintsReal();
        if (subC !is null){
            throw new Exception("constraints in the full system not yet supported",__FILE__,__LINE__);
        }
        return null;
    }
    /+char[] contextId(){ return subContext.contextId; }
    Precision activePrecision() { return subContext.activePrecision; }
    ParticleSys!(Real) refPSysReal() {
        if (_pSysReal is null && subContext.activePrecision==Precision.Real){
            auto nSys=
        }
        return _pSysReal;
    }
    ParticleSys!(LowP) refPSysLowP() { return _pSysLowP; }
    ParticleSys!(Real) pSysReal() { return _pSysReal; }
    ParticleSys!(LowP) pSysLowP() { return _pSysLowP; }
    PSysWriter!(Real) pSysWriterReal(){ return pSysWriter(_pSysReal); }
    PSysWriter!(LowP) pSysWriterLowP() { return pSysWriter(_pSysLowP); }
    void pSysWriterRealSet(PSysWriter!(Real)p) { assert(_pSysReal!is null); _pSysReal[]=p; }
    void pSysWriterLowPSet(PSysWriter!(LowP)p) { assert(_pSysLowP!is null); _pSysLowP[]=p; }
    ConstraintGen constraintGen(){
        auto cReal=constraintsReal();
        if (cReal!is null){
            return cReal.constraintGen();
        }
        auto cLowP=constraintsLowP();
        if (cLowP !is null){
            return cLowP.constraintGen();
        }
        return null;
    }
    
    NotificationCenter nCenter()  { return _nCenter; }
    HistoryManager!(LowP) posHistory() { return _posHistory; }
    ChangeLevel changeLevel() { return _changeLevel; }
    void changeLevelSet(ChangeLevel c) { _changeLevel=c; }
    Real potentialEnergy(){
        mixin(withPSys("return pSys.dynVars.potentialEnergy;"));
    }
    void posSet(SegmentedArray!(Vector!(Real,3)) newPos){
        mixin(withPSys("pSys.dynVars.x.pos[]=newPos;"));
    }
    SegmentedArray!(Vector!(Real,3)) pos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.x.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.x.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.x.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    void dposSet(SegmentedArray!(Vector!(Real,3)) newDpos){
        mixin(withPSys("pSys.dynVars.dx.pos[]=newDpos;"));
    }
    SegmentedArray!(Vector!(Real,3)) dpos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.dx.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.dx.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.dx.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    void mddposSet(SegmentedArray!(Vector!(Real,3)) newMddpos){
        mixin(withPSys("pSys.dynVars.mddx.pos[]=newMddpos;"));
    }
    SegmentedArray!(Vector!(Real,3)) mddpos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.mddx.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.mddx.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.mddx.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    SysStruct sysStruct(){
        SysStruct res;
        mixin(withPSys("res=pSys.sysStruct;"));
        return res;
    }
    void changedDynVars(ChangeLevel changeLevel,Real diff){
        if (changeLevel<_changeLevel) _changeLevel=changeLevel;
        if (diff>maxChange) maxChange=diff;
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
        // maxChange=0.0; changeLevel=ChangeLevel.SmoothPosChange;
    }
    
    void setup(LinearComm pEnv,CharSink log){ }
    
    /// called automatically after creation, but before any energy evaluation
    /// should be called before working again with a deactivated calculator
    void activate(){
        nCenter.notify("willActivateContext",Variant(this));
    }
    /// call this to possibly get rid of all caches (i.e. before a pause in the calculation)
    void deactivate(){
        nCenter.notify("willDeactivateContext",Variant(this));
    }
    /// call this to remove the context (after all calculations with this are finished)
    void giveBack(){
        nCenter.notify("willGiveBackContext",Variant(this));
    }
    /// tries to stop a calculation in progress. Recovery after this is not possible
    /// giveBack should still be called
    void stop(){}
    /// the method of this calculator (this might be different from the method that was used to create this
    /// as that might have been wrapped)
    Method method(){
        assert(0,"to implement in subclasses");
    }
    /// stores the history somewhere and returns an id to possibly recover that history at a later point
    /// this is just an optimization, it does not have to do anything. If implemented then
    /// method.getCalculator, .dropHistory and .clearHistory have to be implemented accordingly
    ubyte[]storeHistory(){ return []; }
    /// exposes (publish/vends) this object to the world
    void publish(){
    }
    mixin(rpcMixin("dchem.CalcContext", "CalculationContext",calcCtxMethodsStr));
    DefaultVendor vendor;
    /// url to access this from other processes (only as CalculationContext)
    char[] exportedUrl(){
        return vendor.proxyObjUrl();
    }
    this(char[] contextId,CharSink log){
        this._contextId=contextId;
        this._logger=log;
        _nCenter=new NotificationCenter();
        // register to the world...
        vendor=new DefaultVendor(this);
        assert(ProtocolHandler.defaultProtocol!is null,"defaultProtocol");
        assert(ProtocolHandler.defaultProtocol.publisher!is null,"publisher");
        ProtocolHandler.defaultProtocol.publisher.publishObject(vendor,"CalcContext"~contextId,true,Publisher.Flags.Public);
    }+/
}


