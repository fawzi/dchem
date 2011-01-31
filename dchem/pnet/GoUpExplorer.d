module dchem.pnet.GoUpExplorer;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;
import dchem.Common;
import blip.container.MinHeap;
import blip.container.HashSet;
import dchem.input.WriteOut;
import blip.math.IEEE;
import blip.sync.Atomic;
import blip.util.LocalMem;
import blip.parallel.mpi.MpiModels;
import dchem.pnet.EmptyExplorer;
import dchem.pnet.PointEvalOp;

class GoUpExplorerDef:ExplorerGen{
    long nEval=long.max;
    this(){
    }
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.MinEExplorer",`
    nEval : number of evaluations to perform (long.max)`));
    mixin printOut!();
    bool verify(CharSink s){
        return true;
    }
    
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos){
        return explorerReal(silos);
    }
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos){
        return explorerLowP(silos);
    }
    ExplorerI!(Real) explorerReal(LocalSilosI!(Real) silos){
        auto res=new GoUpExplorer!(Real)(this,silos);
        return res;
    }
    ExplorerI!(LowP) explorerLowP(LocalSilosI!(LowP) silos){
        auto res=new GoUpExplorer!(LowP)(this,silos);
        return res;
    }
}

// exploration:
// find smallest energy still "free"
// evaluate next direction of it
// apply constraints, if movement is too small declare it as fully visited and neighbor 1
// bcast point as explored
// possibly wait for non collision confirmation
// start calculation
// if too close to existing points stop calculation???
// when calculation is finished
// if mainpoint:
//    gradient -> orient neighbors, compile visited flags, perform topology analysis (attractor, min,max,...)
//    second deriv check
// else store energy, first deriv check??

/// an object that can offer new points to explore
/// actually should not inherit from ExplorationObserverI, but this way we avoid multiple inheritance bugs
class GoUpExplorer(T):EmptyExplorer!(T){
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    HashSet!(Point) removedPoints;
    void delegate(ExplorerI!(T)) available;
    LocalSilosI!(T) silos;
    long leftEval;
    long availEval;
    
    GoUpExplorerDef input;
    this(GoUpExplorerDef input,LocalSilosI!(T)silos){
        this.input=input;
        this.toExploreMore=new MinHeapSync!(PointAndEnergy)();
        this.removedPoints=new HashSet!(Point);
        this.silos=silos;
        this.leftEval=input.nEval;
    }
    
    override string name(){
        return "MinEExplorer_"~input.myFieldName;
    }

    /// adds energy for a local point and bCasts addEnergyEval
    override void addEnergyEvalLocal(SKey key,Point p,Real energy){
        void delegate(ExplorerI!(T)) availableAtt;
        synchronized(this){
            toExploreMore.push(PointAndEnergy(p,energy));
            availableAtt=available;
            available=null;
        }
        if (availableAtt!is null && availEval<leftEval+10){
            availableAtt(this);
        }
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    override void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
        if (silos.hasKey(owner) && !isNaN(pos.potentialEnergy)){
            assert((flags&GFlags.EnergyInfo)==GFlags.EnergyKnown,"non NAN energy, but not EnergyKnown");
            if ((flags&(GFlags.DoNotExplore|GFlags.FullyExplored))==0){
                PointAndEnergy p;
                p.point=point;
                p.energy=pos.potentialEnergy;
                toExploreMore.push(p);
            }
        }
    }
    
    /// finished exploring the given point (remove it from the active points)
    override void finishedExploringPoint(SKey k,Point p,SKey owner){
        if (silos.hasKey(owner))
            rmPoint(p);
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    override void publishCollision(SKey e,Point p){
        rmPoint(p);
    }
    
    // ExplorerI(T)
    
    /// if it had an operation to evaluate, is called on all core silos in parallel
    ReturnFlag nextOp(void delegate(ExplorerI!(T)) available,int req){
        bool availableCalled=false;
        auto nSilos=silos.nextSilosRank(req);
        auto leftEvalNow=atomicAdd(leftEval,-1UL);
        if (leftEvalNow<=0){
            return ReturnFlag.NoOp; // end this explorer
        }
        Point res;
        PointAndEnergy[] pointToDrop;
        scope(exit){
            foreach(p;pointToDrop){
                toExploreMore.push(p);
            }
        }
        while(true) {
            ubyte[512] buf;
            auto lMem=LocalMem(buf);
            PointAndEnergy pe;
            auto paraEnv=silos.paraEnv();
            PointAndEnergy[] lowestPoints=lMem.allocArr!(PointAndEnergy)(paraEnv.dim);
            int iMin;
            for (int iter=0;iter<10;++iter){
                while (toExploreMore.popNext(pe)){
                    synchronized(this){
                        if ((!removedPoints.contains(pe.point)) && 
                            (pointToDrop.length==0 || find(pointToDrop,pe)==pointToDrop.length)) break;
                    }
                }
                if (pe.point.isValid){
                    synchronized(this){
                        toExploreMore.push(pe);
                        this.available=null;
                    }
                }
                if (!availableCalled && pe.point.isValid){
                    if (paraEnv.myRank==0){
                        available(this); // speculatively give back explorer
                    }
                }
                auto dim=paraEnv.dim();
                mpiAllGatherT(paraEnv,pe,lowestPoints);
                if (lowestPoints[0].point.isValid){
                    availableCalled=true;
                }
                double minE=double.min;
                iMin=0;
                foreach(i,pEn;lowestPoints){
                    if (minE>=pEn.energy && pEn.point.isValid){
                        minE=pEn.energy;
                        iMin=i;
                    }
                }
                if (lowestPoints[iMin].point.isValid) break;
            }
            if (!lowestPoints[iMin].point.isValid) {
                synchronized(this){
                    if (!availableCalled){
                        assert(this.available==null || this.available is available,"invalid available value");
                        this.available=available;
                        this.availEval=leftEvalNow;
                    }
                    lMem.deallocArr(lowestPoints);
                    return ReturnFlag.SkipOp; // needs to wait...
                }
            }
            pe=lowestPoints[iMin];
            lMem.deallocArr(lowestPoints);
            if (!availableCalled) {
                availableCalled=true;
                if (paraEnv.myRank==0) available(this); // speculatively call available
            }
            PointAndDir pDir;
            if (iMin==paraEnv.myRank){
                auto lowestPoint=silos.mainPointL(pe.point);
                pDir=lowestPoint.exploreNext();
                if (iMin!=nSilos) {
                    paraEnv[nSilos].send((cast(ubyte*)&pDir)[0..typeof(pDir).sizeof],req);
                }
            }
            EvalOp!(T) newOp;
            if (nSilos==paraEnv.myRank){
                if (iMin!=nSilos) {
                    auto resT=(cast(ubyte*)&pDir)[0..typeof(pDir).sizeof];
                    paraEnv[iMin].recv(resT,req);
                }
                auto lowestPoint=silos.createLocalPoint(pe.point,ev_time());
                scope(exit) { lowestPoint.release(); }
                if (! pDir.point.isValid){ // pe is already fully explored
                    res=Point(1);
                } else if (pDir.dir==invalidDir){
                    res=Point(2);
                } else if (pDir.dir==0){
                    res=pDir.point;
                    newOp=new PointEvalOp!(T)(pDir.point,true);
                    silos.registerPendingOp(newOp);
                    silos.addEvalOp(SKeyVal.Master,newOp,false);
                } else {
                    auto pSysNew=lowestPoint.createPosInDir(pDir.dir);
                    auto newP=silos.newPointAt(pSysNew.dynVars.x,Point(0));
                    if(lowestPoint.acceptNewDirection(newP.point,pDir.dir)){
                        newP=silos.bcastPoint(newP);
                        res=newP.point;
                        newOp=new PointEvalOp!(T)(newP.point);
                        silos.registerPendingOp(newOp);
                        silos.addEvalOp(SKeyVal.Master,newOp,false);
                    } else {
                        newP.drop();
                        res=Point(1);
                    }
                    pSysNew.release();
                }
            }
            mpiBcastT(paraEnv,res,nSilos);
            if (res.isValid){
                return ReturnFlag.LocalOp;
            } else if (res.data==2){
                if (iMin==paraEnv.myRank){
                    pointToDrop~=pe;
                    PointAndEnergy pN;
                    if (toExploreMore.popNext(pN) && pN.point != pe.point)
                        toExploreMore.push(pN);
                }
            } else if (iMin==paraEnv.myRank){
                removedPoints.add(pe.point);
            }
        }
    }
    /// called when an evaluation fails
    override void evaluationFailed(SKey k,Point p){
        rmPoint(p);
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    override bool speculativeGradientLocal(Point p,Real energy){
        synchronized(toExploreMore){
            auto l=toExploreMore.length;
            if (l>0 && toExploreMore.peek.energy>=energy){
                return true;
            }
        }
        return false;
    }
    void rmPoint(Point p){
        removedPoints.add(p);
    }
}
