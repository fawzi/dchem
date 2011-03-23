/// a looper (i.e object that handles the explicit symmetries), that loops on all permutations
/// of groups of particles
module dchem.calculator.SwapSymmNeighLooper;
import dchem.sys.ParticleRange;
import dchem.Common;
import blip.serialization.Serialization;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.sys.SubMapping;
import blip.BasicModels;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.container.BulkArray;
import dchem.input.RootInput;
import dchem.sys.DynVars;
import dchem.calculator.CalculatorModels;
import dchem.calculator.Calculator;
import blip.parallel.mpi.Mpi;
import blip.parallel.smp.Wait;
import dchem.calculator.ProcContext;
import blip.util.LocalMem;


/// uses the SymmNeighLooper of the full system
class SwapSymmNeighLooperGen:SymmNeighLooperGen{
    bool useFirstNeigh=true;
    InputField[] groups;
    mixin(serializeSome("SwapSymmNeighLooper",`handles permutational symmetry of groups of particles`,
    `useFirstNeigh: if first image convention should be used for periodic dimensions (default is true)
    groups: the groups that should be swapped (ParticleRange), all permutations within one group are checked`));
    mixin myFieldMixin!();
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        foreach (g;groups){
            if (g is null){
                dumper(s)("elements of groups should not be null in field ")(myFieldName)("\n");
                res=false;
            } else if ((cast(ParticleRange)g) is null){
                dumper(s)("elements of groups should be ParticleRanges null in field ")(myFieldName)("\n");
                res=false;
            }
        }
        return res;
    }
    SymmNeighLooper symmNeighLooperForContext(CalculationContext ctx){
        return new SwapSymmNeighLooper(this,ctx);
    }
}

/// an object that loops on all structures equivalent to "neigh" up to a permutation in the particles listed in
/// groups that are within epsilon of pSys
///
/// to do: use neigh list (loop only on neighs *and* not yet matched), cache distances
class SwapSymmNeighLooper:SymmNeighLooper{
    SwapSymmNeighLooperGen input;
    CalculationContext ctx;

    this(SwapSymmNeighLooperGen input,CalculationContext ctx){
        this.input=input;
        this.ctx=ctx;
    }
    int loopOnNeighWithinT(T)(ParticleSys!(T)pSys,DistOps distOps,DynPVector!(T,XType)neigh,T epsilon,
        int delegate(ref DynPVector!(T,XType))loopBody)
    {
        // sequential, non cached version
        size_t lTot=0;
        auto groupIdx=new BulkArray!(PIndex)[](input.groups.length);
        auto allIdx=pSys.particlePropertiesPools.particlePropertyT!(int)();
        allIdx[]=0;
        foreach(iGroup,groupF;input.groups){
            size_t lG=0;
            foreach(block;groupF.contentT!(ParticleRange)().loopOn(pSys.sysStruct)){
                lG+=block.length;
                foreach(pIdx;block){
                    (*allIdx.ptrI(pIdx,0))+=1;
                }
            }
            auto gAtt=BulkArray!(PIndex)(lG);
            groupIdx[iGroup]=gAtt;
            size_t ii=0;
            foreach(block;groupF.contentT!(ParticleRange)().loopOn(pSys.sysStruct)){
                foreach(pIdx;block){
                    gAtt[ii++]=pIdx; // randomize???
                }
            }
            lTot+=lG;
        }
        Real distBase=0;
        // do non swapping els
        auto a=pSys.dynVars.x;
        foreach(iRep,pk,lk,v;allIdx.sLoop){
            T[1] dists;
            switch(v){
            case 0:
                distOps.rDistOneToN(pSys,pk.kind,
                    bulkArray(a.pos[pk]),bulkArray(a.orient[pk]),bulkArray(a.dof[pk]),
                    bulkArray(neigh.pos[pk]),bulkArray(neigh.orient[pk]),bulkArray(neigh.dof[pk]),
                    dists);
                distBase+=dists[0];
                break;
            case 1:
                break;
            default:
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("particle ")(pk)(" is present in several permutation groups (")(v)(")");
                }),__FILE__,__LINE__);
            }
            if (distBase>epsilon) return 0;
        }
        Real distAtt=distBase;
        auto configAtt=neigh.dup();
        int res=0;
        // swap dists (sequential, non recursive impl without neighs)
        void swapInGroup(size_t iGroup){
            ubyte[512] buf;
            auto lMem=LocalMem(buf);
            auto fromIdx=groupIdx[iGroup];
            size_t nElGroup=fromIdx.length;
            auto toIdx=fromIdx.dup();
            auto idxAtLevelTo=lMem.allocArr!(size_t)(nElGroup);
            auto distLevel=lMem.allocArr!(Real)(nElGroup);
            scope(exit){
                toIdx.guard.release();
                lMem.deallocArr(idxAtLevelTo);
                lMem.deallocArr(distLevel);
            }
            idxAtLevelTo[0]=0;
            size_t ilevel=0,iTo=0;
            while(res==0){
                // check current dist
                T[1] dists;
                auto pk1=fromIdx[ilevel];
                auto pk2=toIdx[iTo];
                distOps.rDistOneToN(pSys,pk1.kind, // should probably optimize this...
                    bulkArray(a.pos[pk1]),bulkArray(a.orient[pk1]),bulkArray(a.dof[pk1]),
                    bulkArray(neigh.pos[pk2]),bulkArray(neigh.orient[pk2]),bulkArray(neigh.dof[pk2]),
                    dists);
                Real distNew=dists[0];
                if (distNew+distAtt<epsilon){ // go deeper
                    if (iTo!=ilevel){
                        auto tmp=toIdx[iTo];
                        toIdx[iTo]=toIdx[ilevel];
                        toIdx[ilevel]=tmp;
                    }
                    idxAtLevelTo[ilevel]=iTo;
                    distLevel[ilevel]=distAtt;
                    distAtt+=distNew;
                    ++ilevel;
                    if (ilevel!=nElGroup){
                        iTo=ilevel;
                    } else {
                        // do op
                        // mkSwapsIn neigh
                        for (size_t iPerm=0;iPerm<nElGroup;++iPerm){ // could be parallel
                            (configAtt.pos[fromIdx[iPerm]])[]=neigh.pos[toIdx[iPerm]];
                            (configAtt.orient[fromIdx[iPerm]])[]=neigh.orient[toIdx[iPerm]];
                            (configAtt.dof[fromIdx[iPerm]])[]=neigh.dof[toIdx[iPerm]];
                        }
                        if (iGroup<groupIdx.length){
                            swapInGroup(iGroup+1);
                        } else {
                            auto rAtt=loopBody(configAtt);
                            if (rAtt!=0) res=rAtt;
                        }
                    }
                } else { // same level
                    while(true){
                        ++iTo;
                        if (iTo==nElGroup){ // backtrack one level
                            if (ilevel==0) return; // end
                            --ilevel;
                            iTo=idxAtLevelTo[ilevel];
                            distAtt=distLevel[ilevel];
                            if (iTo!=ilevel){
                                auto tmp=toIdx[iTo];
                                toIdx[iTo]=toIdx[ilevel];
                                toIdx[ilevel]=tmp;
                            }
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        if (groupIdx.length>0){
            swapInGroup(0);
        } else {
            res=loopBody(neigh);
        }
        configAtt.giveBack();
        return res;
    }
    // aliases don't work reliably
    int loopOnNeighWithin(ParticleSys!(Real)pSys,DistOps distOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody)
    {
        return loopOnNeighWithinT!(Real)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithin(ParticleSys!(LowP)pSys,DistOps distOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody)
    {
        return loopOnNeighWithinT!(LowP)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithinReal(ParticleSys!(Real)pSys,DistOps distOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody)
    {
        return loopOnNeighWithinT!(Real)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithinLowP(ParticleSys!(LowP)pSys,DistOps distOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody)
    {
        return loopOnNeighWithinT!(LowP)(pSys,distOps,neigh,epsilon,loopBody);
    }
}
