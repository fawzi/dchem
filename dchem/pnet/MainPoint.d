/// this module represents 
module dchem.pnet.MainPoint;
import dchem.Common;
import dchem.sys.ParticleSys;
import dchem.sys.DynVars;
import dchem.calculator.CalculatorModels;
import dchem.calculator.Calculator;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import blip.core.Boxer;
import blip.sync.Atomic;
import blip.math.Math:max,abs,pow2,sqrt;
import tango.util.container.more.Heap;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import tango.util.container.HashMap;
import dchem.input.RootInput;
import dchem.input.WriteOut;
import blip.container.Deque;
import dchem.pnet.PNetModels;
import dchem.pnet.DirArray;
import blip.container.Pool;
import blip.container.Cache;
import blip.core.Array;
import blip.core.Traits:isNullT;
import blip.parallel.smp.WorkManager;
import dchem.util.Rotate;
import blip.time.Clock;
import blip.time.Time;
import blip.narray.NArray;
import blip.math.IEEE;
import blip.util.LocalMem;
import blip.io.EventWatcher;

/// structure that collects togheter all points of the same owner (useful to bcast chuncks)
class PointToOwner{
    /// points sorted by owner & point
    Point[] points;
    /// index of the point in the original array
    size_t[] idx;
    /// list of owners
    SKey[] owner;
    /// start of the points of each owner (+ end, owner[i] has points[starts[i]..starts[i+1]])
    size_t[] starts;

    mixin(descSome("PointToOwner",`points|idx|owner|starts`));
    
    version(TrackPointToOwner){
        final void logMsg(void delegate(CharSink) logger){
            sinkTogether(sout.call,delegate void(CharSink s){
                dumper(s)("<PointToOwner_")(cast(void*)this)(">");
                logger(s);
                dumper(s)("</PointToOwner_")(cast(void*)this)(">");
            });
        }
    } else {
        final void logMsg(void delegate(CharSink) logger){}
    }
    /// internal help structure (to reorder stuff)
    struct PIdxOwn{
        Point point;
        size_t idx;
        SKey owner;
        mixin(serializeSome("dchem.PointToOwner.PIdxOwn",``,`point|idx|owner`));
        mixin printOut!();
    }
    /// returns a PointToOwner that contains the points extracted with getPnt ordered by owner and point
    /// if compressPts is true (default) then duplicate points are removed
    /// if keepIdx is true (the default) the indexes of the points in the origina array are kept
    this(SKey delegate(Point) ownerMap,Point delegate(size_t i) getPnt, size_t nPts, bool compressPts=true,bool keepIdx=true){
        if (nPts==0) return;

        PIdxOwn[128] buf;
        auto pts=lGrowableArray(buf,0);
        scope(exit){ pts.deallocData(); }
        size_t i=0,j=1;
        {
            PIdxOwn pNew;
            pNew.point=getPnt(0);
            pNew.idx=0;
            pNew.owner=ownerMap(getPnt(0));
            pts(pNew);
        }
        while(j<nPts){
            if ((!compressPts) || pts[i].point!=getPnt(j)){
                PIdxOwn pNew;
                pNew.point=getPnt(j);
                pNew.idx=j;
                pNew.owner=ownerMap(getPnt(j));
                ++i;
                pts(pNew);
            }
            ++j;
        }
        ++i;
        auto points=pts.data;
        version(TrackPointToOwner){
            logMsg((CharSink s){
                dumper(s)("nonOrdPoints:")(points);
            });
        }
        sort(points,delegate bool(PIdxOwn a,PIdxOwn b){ return a.owner<b.owner || (a.owner==b.owner && a.point<b.point); });
        version(TrackPointToOwner){
            logMsg((CharSink s){
                dumper(s)("ordPoints:")(points);
            });
        }
        if (compressPts){
            points=compress(points,delegate bool(PIdxOwn a,PIdxOwn b){ return a.point==b.point; });
            version(TrackPointToOwner){
                logMsg((CharSink s){
                    dumper(s)("compressed Points:")(points);
                });
            }
        }
        this.points=new Point[](points.length);
        for (size_t ii=0;ii<points.length;++ii){
            this.points[ii]=points[ii].point;
        }
        this.idx=new size_t[](points.length);
        for (size_t ii=0;ii<points.length;++ii){
            this.idx[ii]=points[ii].idx;
        }
        size_t nOwners=1;
        for (size_t ii=1;ii<points.length;++ii){
            if (points[ii-1].owner!=points[ii].owner)++nOwners;
        }
        this.starts=new size_t[](nOwners+1);
        this.owner=new SKey[](nOwners);
        this.owner[0]=points[0].owner;
        this.starts[0]=0;
        j=0;
        for (size_t ii=1;ii<points.length;++ii){
            if (points[ii-1].owner!=points[ii].owner){
                this.owner[++j]=points[ii].owner;
                this.starts[j]=ii;
            }
        }
        ++j;
        assert(j==nOwners);
        this.starts[j]=points.length;
    }
    /// number of owners
    size_t nOwners(){
        return owner.length;
    }
    /// gets the i-th owner and its points
    bool ptsI(size_t i,out SKey owner,out Point[]points){
        if (i<this.owner.length){
            owner=this.owner[i];
            points=this.points[starts[i]..starts[i+1]];
            return true;
        }
        return false;
    }
    /// gets the i-th owner and its points and indexes
    bool ptsI(size_t i,out SKey owner,out Point[]points,out size_t[] idx){
        if (i<this.owner.length){
            owner=this.owner[i];
            points=this.points[starts[i]..starts[i+1]];
            if (this.idx.length!=0)
                idx=this.idx[starts[i]..starts[i+1]];
            return true;
        }
        return false;
    }
    /// gives back the memory used by this structure
    void deallocData(){
        if (points !is null) delete points;
        if (idx!is null) delete idx;
        if (starts!is null) delete starts;
        if (owner!is null) delete owner;
    }
    /// helper structure for iterations
    struct Iter{
        PointToOwner pOwn;
        size_t iOwn;
        Iter* next;
        PoolI!(Iter*) pool;
        void delegate(Iter*) action;
        
        SKey owner(){
            return pOwn.owner[iOwn];
        }
        Point[] points(){
            return pOwn.points[pOwn.starts[iOwn]..pOwn.starts[iOwn+1]];
        }
        size_t[] idx(){
            return pOwn.idx[pOwn.starts[iOwn]..pOwn.starts[iOwn+1]];
        }
        size_t localLb(){
            return pOwn.starts[iOwn];
        }
        size_t localUb(){
            return pOwn.starts[iOwn+1];
        }
        void giveBack(){
            pool.giveBack(this);
        }
        void doAction(){
            action(this);
        }
    }
    /// helper structure for parallel loops
    struct PLoop{
        PointToOwner pOwn;
        int opApply(int delegate(ref Iter) loopBody){
            auto pool=cachedPoolNext(function Iter*(PoolI!(Iter*)p){
                auto res=new Iter;
                res.pool=p;
                return res;
            });
            int res;
            Exception exception;
            void doIter(Iter *i){
                if (res!=0 || exception!is null) return;
                try{
                    auto rAtt=loopBody(*i);
                    if (rAtt!=0) res=rAtt;
                } catch (Exception e){
                    exception=e;
                }
            }
            Task("PointToOwnerPLoop",delegate void(){
                auto nOwn=pOwn.owner.length;
                for (size_t iOwn=0;iOwn<nOwn;++iOwn){
                    auto nIter=pool.getObj();
                    nIter.pOwn=pOwn;
                    nIter.action=&doIter;
                    nIter.iOwn=iOwn;
                    Task("PointToOwnerPLoopIter",&nIter.doAction).appendOnFinish(&nIter.giveBack).autorelease.submitYield();
                }
            }).autorelease.appendOnFinish(&pool.rmUser).executeNow();
            if (exception!is null) throw new Exception("exception in PointToOwnerPLoop",__FILE__,__LINE__,exception);
            return res;
        }
    }
    /// sequential loop
    int opApply(int delegate(ref Iter) loopBody){
        Iter nIter;
        for (size_t iOwn=0;iOwn<owner.length;++iOwn){
            nIter.iOwn=iOwn;
            auto rAtt=loopBody(nIter);
            if (rAtt!=0) return rAtt;
        }
    }
    /// parallel loop
    PLoop pLoop(){
        PLoop res;
        res.pOwn=this;
        return res;
    }
}

/// structure to describe a GFlags change of a point
struct GFlagsChange{
    uint oldGFlags;
    uint newGFlags; // avoid storing??
    Point point;
    long refId;
    mixin(serializeSome("dchem.GFlagsChange",`describes a change of flags of a point`,"oldGFlags|newGFlags|point|refId"));
    mixin printOut!();
}

/// structure to describe new neighbors
struct PointNeighbors(T){
    Point point;
    PointAndDir[] neighbors;
    DirDistances!(T)[] dirDist;
    Time time;
    mixin(serializeSome("dchem.PointNeighbors!("~T.stringof~")","describes a new neighbor","point|neighbors|dirDist"));
    mixin printOut!();
}

/// a main evaluation point
class MainPoint(T):MainPointI!(T){
    /// the local context of this point (it is not the main owner if this is a LocalCopy)
    LocalSilosI!(T) localContext(){ return _localContext; }
    LocalSilosI!(T) _localContext;
    /// identification of this point
    Point point() { return _point; }
    Point _point;
    /// attractor of this point (attractors are identified by their minima)
    Attractor attractor(){ return _attractor; }
    void attractor(Attractor newA){
        bool minChanged=false;
        synchronized(this){
            if ((isNaN(_attractor.energyThroughPoint) || newA.energyThroughPoint<=_attractor.energyThroughPoint) &&
                (_attractor.throughPoint!=newA.throughPoint || _attractor.idThroughPoint<newA.idThroughPoint)){
                _attractor.energyThroughPoint=newA.energyThroughPoint;
                ++(_attractor.id);
            }
        }
    }
    Attractor _attractor;
    
    /// position of the point (energy, derivatives,... might be invalid, only real positions have to be valid)
    ParticleSys!(T) pos(){ return _pos; }
    ParticleSys!(T) _pos;
    /// direction of the minimum in the dual space with norm 1 wrt. euclidean norm
    /// (frame of reference for minimization), valid only if this point is a starting point for further exploration
    DynPVector!(T,DualDxType) minDir(){ return _minDir; }
    DynPVector!(T,DualDxType) _minDir;
    /// bit array of the directions that have been explored (is of length null all direction have been explored)
    /// this stores only real exploration directions, directions that a given method blocks should be stored in
    /// a separated DirArray
    FlagsArray exploredDirs(){ return _exploredDirs; }
    FlagsArray _exploredDirs;
    /// neighbors, and the direction with respect to this point in which they are
    GrowableArray!(PointAndDir) neighbors(){ return &_neighbors; }
    LocalGrowableArray!(PointAndDir) _neighbors;
    /// neighbor distances, stored only if requested. 6 numbers for each neighbor:
    /// dualDist, cartDist, cosDual, rDistDual, cartCos, cartRDist
    GrowableArray!(DirDistances!(T)) neighDistances(){ return &_neighDistances; }
    LocalGrowableArray!(DirDistances!(T)) _neighDistances;
    /// scale of mindir to recover the dual gradient (useful??)
    T minDirScale(){ return _minDirScale; }
    T _minDirScale;
    /// exploration size for this point (used to establish neighbors,...) 
    T explorationSize(){ return _explorationSize; }
    T _explorationSize;
    /// bit-or of GFlags of the current point
    uint gFlags(){ return _gFlags; }
    uint _gFlags;
    size_t refCount=1;
    PoolI!(MainPoint) pool;
    
    void logMsg(void delegate(CharSink)logger){
        localContext.logMsg(delegate void(CharSink sink){
            auto s=dumper(sink);
            s("<MainPointLog");
            s(" point=")(point.data);
            s(" addr=@")(cast(void*)this);
            if (isLocalCopy) s(" copy=1");
            s(">\n  ");
            indentWriter("  ",sink,logger);
            s("\n</MainPointLog>");
        });
    }
    SKey owner(){
        return localContext.ownerOfPoint(point);
    }
    /// if this point is a local copy, and not the "main" point
    bool isLocalCopy(){
        return (gFlags&GFlags.LocalCopy)!=0;
    }
    /// if the point is explorable
    bool isExplorable(){
        auto flags=gFlags;
        return (flags&GFlags.EnergyInfo)==GFlags.EnergyEvaluated &&
            (flags&(GFlags.DoNotExplore|GFlags.FullyExplored|GFlags.FullyEvaluated|GFlags.OldApprox))==0;
    }
    bool isValid(){
        return (gFlags&GFlags.OldApprox)==0;
    }
    bool isPublic(){
        return (gFlags&GFlags.PointBcastStatus)!=0;
    }
    /// if the frame of reference is set
    bool hasFrameOfRef(){
        return (gFlags&GFlags.HasRefFrame)!=0;
    }
    void bcastLevel(int i){
        uint oldGFlags,newGFlags;
        synchronized(this){
            oldGFlags=_gFlags;
            switch(i){
            case 0:
                _gFlags|=GFlags.PointBcasted;
                break;
            case 1:
                _gFlags|=GFlags.PointNeighBcasted;
                break;
            default:
                assert(0,"unexpected bcastLevel");
            }
            newGFlags=_gFlags;
        }
        if (oldGFlags!=newGFlags) notifyGFlagChange(oldGFlags);
    }
    
    // *** topology, utility methods
    /// probability of having a minimum close by (utility method)
    Prob minimum(){ return minimumForGFlags(gFlags); }
    /// probability of having a critical point close by (utility method)
    Prob criticalPoint(){ return criticalPointForGFlags(gFlags); }
    /// probability of having a saddle point (utility method)
    Prob saddlePoint(){ return saddlePointForGFlags(gFlags); }
    /// probability that close to this point there is a special (critical) point (utility method)
    Prob specialPoint(){ return specialPointForGFlags(gFlags); }
    /// returns the main that this point can be considered as (utility method)
    PointType pointType(){ return pointTypeForGFlags(gFlags); }
    
    /// this point as dried point
    DriedPoint!(T) driedPoint(){
        DriedPoint!(T) res;
        res.point=point;
        res.pos=pSysWriter(pos);
        res.minDir=dynPVectorWriter(minDir);
        res.exploredDirs=exploredDirs;
        res.neighbors=_neighbors;
        res.attractor=attractor;
        res.minDirScale=minDirScale;
        res.explorationSize=explorationSize;
        res.gFlags=_gFlags;
        return res;
    }
    /// if given the neighbor are taken into ownership, and freed by this object.
    this(LocalSilosI!(T) localContext, PoolI!(MainPoint) p,PointAndDir[] neighbors=[]){
        this.refCount=1;
        this._localContext=localContext;
        this.pool=p;
        this._pos=localContext.refPos.dup(PSDupLevel.EmptyDyn);
        this._gFlags=GFlags.None;
        this._explorationSize=cast(T)localContext.explorationSize();
        this._minDirScale=T.init;
        //this._minDir.clear();
        this._exploredDirs=new FlagsArray(ndirs);
        this._neighbors=lGrowableArray(neighbors,neighbors.length,GASharing.Global);
    }
    
    void clear(){
        if (pos!is null) this.pos.dynVars.deallocData();
        this._gFlags=GFlags.None;
        this._explorationSize=cast(T)localContext.explorationSize();
        this._minDirScale=T.init;
        this._minDir.clear();
        this._exploredDirs.clearData();
        this._neighbors.clearData();
        this._neighDistances.clearData();
    }
    /// constructor
    this(LocalSilosI!(T) localContext,Point point,ParticleSys!(T) pos,uint gFlags,
        T explorationSize,DynPVector!(T,DualDxType) minDir,T minDirScale=T.init,
        FlagsArray exploredDirs=null,PointAndDir[] neighbors=null,PoolI!(MainPoint)pool=null)
    {
        this.refCount=1;
        this._localContext=localContext;
        this._point=point;
        this._pos=pos;
        if (this.pos is null){
            this._pos=localContext.refPos.dup(PSDupLevel.EmptyDyn);
        }
        this._gFlags=gFlags;
        this._explorationSize=explorationSize;
        if (this.explorationSize<=0) this._explorationSize=localContext.explorationSize;
        this._minDir=minDir;
        this._minDirScale=minDirScale;
        this._exploredDirs=exploredDirs;
        this.pool=pool;
        if (this._exploredDirs is null){
            this._exploredDirs=new FlagsArray(ndirs);
        }
        this._neighbors=lGrowableArray(neighbors,neighbors.length,GASharing.Global);
    }
    this(LocalSilosI!(T) localContext,Point point,ParticleSys!(T) pos,uint gFlags=GFlags.None,
        T explorationSize=-1)
    {
        DynPVector!(T,DualDxType) minD;
        this(localContext,point,pos,gFlags,explorationSize,minD);
    }
    /// atomic cas on the flags of this point
    uint gFlagsAtomicCAS(uint newVal,uint oldVal){
        uint oVal;
        synchronized(this){
            oVal=_gFlags;
            if (oVal==oldVal){
                _gFlags=newVal;
            }
        }
        if (oVal==oldVal&&newVal!=oldVal) notifyGFlagChange(oVal);
        return oVal;
    }
    /// atomic op on the flags of this point
    uint gFlagsAtomicOp(uint delegate(uint) op){
        uint oldV;
        synchronized(this){
            oldV=_gFlags;
            _gFlags=op(oldV);
        }
        if (_gFlags!=oldV) notifyGFlagChange(oldV);
        return oldV;
    }
    /// explores the next direction, immediately marks it as in exploration
    /// returns direction 0 if the gradient has to be evaluated, an invalid direction if the current point is in evaluation,
    /// and an invalid point only if all direction are explored/ the point should not be explored
    /// if lastIsLast is true then the last direction (gradient) is explored as last
    PointAndDir exploreNext(FlagsArray methodDirs=null,bool lastIsLast=true){
        auto exploreFlags=DirFlags.Explored;
        if ((gFlags&GFlags.HasRefFrame)==0) {
            synchronized(this){
                if ((gFlags&(GFlags.InProgress|GFlags.GradientEvaluated))==0){
                    if (exploredDirs.atomicCAS(0,exploreFlags,DirFlags.Free)==DirFlags.Free){
                        _gFlags|=GFlags.GradientInProgress;
                        version(TrackPNet) logMsg(delegate void(CharSink s){
                            dumper(s)("exploreNext will return self Eval");
                        });
                        return PointAndDir(point,0);
                    } else if ((gFlags&GFlags.DoNotExplore)==0){
                        throw new Exception("direction 0 taken, but no gradient eval in progress",__FILE__,__LINE__);
                    }
                }
            }
        }
        if ((gFlags&GFlags.DoNotExplore)!=0){
            version(TrackPNet) logMsg(delegate void(CharSink s){
                dumper(s)("GFlags.DoNotExplore: exploreNext will return 0,invalidDir");
            });
            return PointAndDir(Point(0),invalidDir);
        }
        assert((gFlags&GFlags.PointBcasted)!=0);
        if ((gFlags&GFlags.InProgress)!=0){
            version(TrackPNet) logMsg(delegate void(CharSink s){
                dumper(s)("GFlags.InProgress: exploreNext asks to wait");
            });
            return PointAndDir(this.point,invalidDir); // needs to wait evaluation
        }
        if ((gFlags&GFlags.PointNeighBcasted)==0){
            version(TrackPNet) logMsg(delegate void(CharSink s){
                dumper(s)("Not yet fully bcasted neighs: exploreNext asks to wait");
            });
            return PointAndDir(this.point,invalidDir); // needs to wait evaluation
        }
        if ((gFlags&GFlags.FullyExplored)==0){
            uint oldGFlags,newGFlags;
            synchronized(this){
                oldGFlags=gFlags;
                if ((gFlags&GFlags.FullyExplored)==0){
                    size_t nextDir=exploredDirs.length;
                    if(lastIsLast && (methodDirs is null || methodDirs.atomicCAS(1,exploreFlags,0)==0) && 
                        exploredDirs.atomicCAS(1,exploreFlags,0)==0) // explore down first
                    {
                        nextDir=1;
                    } else {
                        if (lastIsLast) --nextDir;
                        auto start=localContext.rand.uniformR(nextDir);
                        nextDir=exploredDirs.findFreeAndSet(start,exploreFlags,lastIsLast,methodDirs);
                    }
                    version(TrackPNet) logMsg(delegate void(CharSink s){
                        dumper(s)("exploreNext nextDir:")(nextDir)(" of ")(exploredDirs.length);
                    });
                    if (nextDir!=exploredDirs.length) return PointAndDir(this.point,nextDir);
                    _gFlags|=GFlags.FullyExplored;
                }
                newGFlags=newGFlags;
            }
            if (oldGFlags!=newGFlags) notifyGFlagChange(oldGFlags);
        }
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("exploreNext on fully explored point");
        });
        return PointAndDir(Point(0),0);
    }
    /// returns the number of dimensions
    uint ndim(){
        auto g=pos.dynVars.dVarStruct.dualDxGroup;
        return 3*g.posStruct.length+4*g.orientStruct.length+g.dofStruct.length;
    }
    /// returns the number of directions (counting also 0, the core dir)
    uint ndirs(){
        return 2*ndim+1;
    }
    /// returns the dir value for the given dimension and sign
    uint toDir(uint idim,bool neg){
        assert(idim<ndim,"dim out of bounds");
        if (idim==0 && neg) return ndirs-1;
        return 2*idim+1-(neg?1:0);
    }
    /// transforms a non core (i.e. non 0) dir to dimension and sign
    void fromDir(uint dir,out uint dim,out bool neg){
        assert(dir>0,"special core direction not mappable");
        assert(dir<ndirs,"dir is too large");
        dim=dir/2;
        if(dir==ndirs-1) dim=0;
        neg=((dir&1)==0);
    }
    /// a point in the given direction in the dual space
    DynPVector!(T,DualDxType) createDualDir(uint dir){
        DynPVector!(T,DualDxType) v0;
        if (dir==0) return v0; // return a null vector?? (move two lines down)
        v0=minDir.emptyCopy;
        v0[]=0;
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        T val=(neg?-1:1);
        v0[rDir]=val;
        auto newDir=rotateEiV(0,minDir,v0);
        return newDir;
    }
    /// energy of the current point
    Real energy(){
        synchronized(this){
            return pos.dynVars.potentialEnergy;
        }
    }
    /// energy of the current point
    Real energyError(){
        synchronized(this){
            return pos.dynVars.potentialEnergyError;
        }
    }
    PointEMin pointEMin(){
        PointEMin res;
        synchronized(this){
            res.energy=pos.dynVars.potentialEnergy;
            res.point=point;
            res.minimum=_attractor.minimum;
            res.id=_attractor.id;
        }
        return res;
    }
    /// projects into the correct movement space (back & forth from the direct space).
    /// and gets rid also of the elements in the constraint space
    /// useful to get rid of the components in the null and constraint space
    T projectInDualTSpace(ParticleSys!(T) pSys,DynPVector!(T,DualDxType)dualDeriv){
        T res;
        auto overlap=pSys.maybeDerivOverlap();
        if (overlap!is null){
            scope newDirDirect=pSys.dynVars.dVarStruct.emptyDx();
            pSys.fromDualTSpace(dualDeriv,newDirDirect);
            localContext.constraints().applyDR(pSys,newDirDirect);
            pSys.toDualTSpace(newDirDirect,dualDeriv);
            res=pSys.dotInTSpace(newDirDirect,dualDeriv);
            newDirDirect.giveBack();
        } else {
            localContext.constraints().applyDR(pSys,dualDeriv.toGroup!(DxType)());
            res=dualDeriv.norm2();
        }
        return res;
    }
    
    /// calculates the position exploring from here in the given direction
    /// returns null if the position is *really* too close
    ParticleSys!(T) createPosInDir(uint dir){
        version (TrackPNet) logMsg(delegate(CharSink s){
            dumper(s)("called createPosInDir(")(dir)(")\n");
            dumper(s)("minDir:")(dynPVectorWriter(minDir));
        });
        auto newDir=createDualDir(dir);
        version (TrackPNet) logMsg(delegate(CharSink s){
            dumper(s)("pippo newDir:")(dynPVectorWriter(newDir));
        });
        auto cartesianNorm=projectInDualTSpace(pos,newDir);
        version (TrackPNet) logMsg(delegate(CharSink s){
            dumper(s)("pippo proj newDir:")(dynPVectorWriter(newDir))(" cartesianNorm:")(cartesianNorm);
        });
        auto projNorm=newDir.norm2();
        version (TrackPNet) logMsg(delegate(CharSink s){
            dumper(s)("pippo projNorm:")(projNorm);
        });
        if (cartesianNorm<=localContext.minRealNormSelf0){
            logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" discarded because cartesian norm is too small. Projected norm was ")(projNorm)(".");
            });
            return null; // continue all the same or return pos???
        }
        assert(projNorm>0,"projNorm should not be 0 without cartesian norm being 0 too");
        auto scaleAtt=explorationSize/projNorm; // always normalizing in the non null directions
        /// immediately check expected cartesian norm change and possibly increase length
        /// this could be wrong for periodic directions, but is most likely better than not doing it
        if (cartesianNorm*scaleAtt<localContext.minRealNormSelf1){
            logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" has small cartesian norm, increasing search length");
            });
            scaleAtt=localContext.minRealNormSelf1/cartesianNorm;
        }
        newDir*=scaleAtt;
        
        auto resPSys=pos.dup(PSDupLevel.DynProperties|PSDupLevel.DynPNullDx|PSDupLevel.DynPNullMddx|PSDupLevel.HiddenVars);
        pos.addFromDualTSpace!(T)(newDir,resPSys.dynVars.x);
        newDir.giveBack();
        auto err=localContext.constraints.applyR(resPSys);
        if (err>0.1){
            logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)("discarded because we might have a constraint clash");
            });
            return null;
        }
        resPSys.updateHVars();
        return resPSys;
    }
    
    /// returns if the given direction is acceptable. If it is accepted sets the local directions covered
    /// this is most likely called on a copy of the point, and where newPoint is local...
    bool acceptNewDirection(Point newPoint,uint dir){
        // Point
        if (dir==0){
            throw new Exception("should not call acceptNewDirection when exploring the point itself (direction 0)",
                __FILE__,__LINE__);
        }
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        auto localPoint=localContext.createLocalPoint(newPoint,ev_time());
        scope(exit) { localPoint.release(); }
        auto newPos=localPoint.pos.dynVars.x;
        
        /// verify pos:
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);
        DirDistances!(T)[128] bufDirDist;
        auto dirDist=lGrowableArray(bufDirDist,0);
        scope(exit){
            neighAtt.deallocData();
            dirDist.deallocData();
        }
        bool hadGrad=hasFrameOfRef;
        auto mainDists=addDirsOf(newPos,newPoint,localPoint.explorationSize,neighAtt,dirDist);

        if (mainDists.xDist<localContext.minMoveInternal){
            logMsg(delegate void(CharSink s){
                dumper(s)("norm in the internal coordinates of exploration direction ")(dir)(" is too small in internal coordinates ")(mainDists.xDist)(", discarding evaluation");
            });
            return false;
        }
        // original norm in dual space
        auto origNorm=mainDists.dualDist;
        if (origNorm<=localContext.minNormDualSelf()*explorationSize){
            logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)
                    (" is too small:")(origNorm)(", discarding evaluation");
            });
            return false;
        }
        if (origNorm>localContext.maxNormDual*explorationSize || origNorm<localContext.minNormDual*explorationSize){
            logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)
                    (" is outside the exploration zone with length ")(origNorm/explorationSize)("*")(explorationSize)
                    (", continuing evaluation");
            });
        }
        // norm in cartesian units
        auto normCartesian=mainDists.cartesianDist;
        if (normCartesian<localContext.minRealNormSelf2){
            logMsg(delegate void(CharSink s){
                dumper(s)("norm in real space of exploration in direction ")(dir)
                    (" is too small:")(normCartesian)(", discarding evaluation");
            });
        }
        
        if (findFirstPred(neighAtt.data,delegate bool(PointAndDir p){ return p.dir==dir; })==neighAtt.length){
            auto diff=newPos.dup();
            diff.opBypax(pos.dynVars.x,-1,1);
            auto deriv1=pos.dynVars.dVarStruct.emptyDx();
            deriv1[]=0;
            pos.addToTSpace!(T)(diff,deriv1);
            diff.giveBack();
            // make an optimized version when the overlap is the identity??
            auto deriv1Dual=pos.dynVars.dVarStruct.emptyDualDx();
            pos.toDualTSpace!(T)(deriv1,deriv1Dual);
            scope(exit){
                deriv1Dual.giveBack();
                deriv1.giveBack();
            }
            // go in the minDir frame of reference
            auto deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
            auto deriv2=rotateVEi(minDir.toGroup!(DxType)(),0,deriv1);
            
            // not in the expected direction
            if (neighAtt.length==0){ // not a neighbor...
                logMsg(delegate void(CharSink s){
                    dumper(s)("exploration in direction ")(dir)
                        (" is not a neighbor, continuing evaluation..."); // stop instead?
                });
                return true;
            } else {
                // check the direction
                
                /// checking which direction of the dual space we are really exploring
                volatile T maxVal=0;
                size_t iMax;
                foreach (i,v;deriv2Dual.sLoop){
                    auto vp=abs(v);
                    if (vp>maxVal){
                        maxVal=vp;
                        iMax=i;
                    }
                }
                
                uint dirMax=toDir(iMax,deriv2Dual[iMax]<0);
                bool nonVisited=false;
                for (size_t idir=0;idir<neighAtt.length;++idir){
                    if (neighAtt[idir].dir==0){
                        if (dirMax==dir){
                            logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)
                                    (" is  mostly in the expected direction in the dual space, but ")
                                    (iMax)(", but not in the expected region, continuing evaluation");
                            }); // change, do not visit or require at least a minimum cartesian norm???
                            nonVisited=true;
                        } else {
                            logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)
                                    (" is  very close to the starting point, and not in the expected direction");
                            });
                        }
                    } else {
                        auto actualDirFlags=exploredDirs.atomicCAS(dirMax,DirFlags.Explored,DirFlags.Free);
                        if (actualDirFlags == DirFlags.Free){
                            logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)
                                    (" is not mostly in the expected direction in the dual space, but in dir ")
                                    (neighAtt[idir].dir)(", continuing evaluation declaring also it as visited");
                            });
                            nonVisited=true;
                        } else {
                            logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)
                                    (" is not mostly in the expected direction in the dual space, but in dir ")
                                    (neighAtt[idir].dir)(", and that was already visited");
                            });
                        }
                    }
                }
                if (!nonVisited) {
                    logMsg(delegate void(CharSink s){
                        dumper(s)("exploration in direction ")(dir)
                            (" is not mostly in the expected direction in the dual space, but in directions that ")
                            ("where already visited, discarding it");
                    });
                    return false;
                }
            }
        }
        assert(neighAtt.length!=0);
        addNeighbors(neighAtt.data,dirDist.data,hadGrad);
        return true;
    }
    
    /// internal method add the pre-processed points&dirs to the list of neighbors
    /// and notify neighbors
    bool addNeighbors(PointAndDir[] neighs,DirDistances!(T)[] dirDist,bool hadGradient)
    in {
        auto newPoint=neighs[0].point;
        foreach(p;neighs){
            assert(p.point==newPoint,"expected addition of a single point"); // lift constraint??
        }
    } body {
        version(TrackPNet){
            logMsg(delegate void(CharSink s){
                dumper(s)("addNeighbors(")(neighs)(",")(dirDist)(",")(hadGradient)(")");
            });
        }
        if (neighs.length==0) return false; // return true?
        if (isLocalCopy){
            bool added=localContext.addNeighDirsToLocalPoint(localContext.ownerOfPoint(point),point,neighs,dirDist,hadGradient);
            return added;
        }
        bool gradChanged=false;
        auto newPoint=neighs[0].point;
        synchronized(this){
            /// find if we already added this point...
            auto pos=findFirstPred(neighbors.data,delegate bool(PointAndDir pD){ return pD.point==newPoint; });
            if (pos!=neighbors.data.length){
                return false;
            }
            if ((!hadGradient) && hasFrameOfRef){
                gradChanged=true;
            } else {
                auto nd=ndirs;
                foreach(n;neighs){
                    if (0<n.dir && n.dir<nd){
                        exploredDirs.atomicCAS(n.dir,DirFlags.Explored,DirFlags.Free);
                    }
                }
                neighbors.appendArr(neighs); // duplicates are stored contiguously
                _neighDistances.appendArr(dirDist);
            }
        }
        if (gradChanged){
            version(TrackPNet) logMsg(delegate void(CharSink s){ s("gradChanged"); });
            PointAndDir[128] buf;
            auto neighAtt=lGrowableArray(buf,0);
            DirDistances!(T)[128] bufDirDist;
            auto dirDist1=lGrowableArray(bufDirDist,0);
            scope(exit){
                neighAtt.deallocData();
                dirDist1.deallocData();
            }
            auto localPoint=localContext.createLocalPoint(newPoint,ev_time());
            scope(exit){ localPoint.release(); }
            auto newPos=localPoint.pos.dynVars.x;
            auto dDir=addDirsOf(newPos,newPoint,localPoint.explorationSize,neighAtt,dirDist1);
            return addNeighbors(neighAtt.data,dirDist1.data,true);
        }
        // communicate to neighbors
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("communicate to neighbors the addition of ")(newPoint.data);
        });
        auto neigToNotify=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return neighs[i].point; },neighs.length);
        foreach (iPts;neigToNotify.pLoop){
            localContext.addPointToLocalNeighs(iPts.owner,point,iPts.points);
        }
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("notify PointHasNewNeighbors ")(newPoint.data);
        });
        // notify
        PointNeighbors!(T) pn;
        pn.point=point;
        pn.neighbors=neighs;
        pn.dirDist=dirDist;
        pn.time=Clock.now;
        localContext.nCenter.notify("PointHasNewNeighbors",box(pn));
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("completed addition of ")(newPoint.data);
        });
        return true;
    }

    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point.
    bool checkIfNeighbor(DynPVector!(T,XType)newPos,Point newPoint,T pSize) {
        // ignore self
        if (newPoint==point) return false; // should return true???
        // ignore if the point is already in the neighbors
        synchronized(this){
            if (findFirstPred(neighbors.data,delegate bool(PointAndDir p){ return p.point==newPoint; })!=neighbors.length){
                /// already added
                return true;
            }
        }
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);
        DirDistances!(T)[128] bufDirDist;
        auto dirDist=lGrowableArray(bufDirDist,0);
        scope(exit){
            neighAtt.deallocData();
            dirDist.deallocData();
        }
        bool hasGrad=hasFrameOfRef;
        auto res=addDirsOf(newPos,newPoint,pSize,neighAtt,dirDist);
        if (neighAtt.length>0){
            addNeighbors(neighAtt.data,dirDist.data,hasGrad);
            return true;
        }
        return false;
    }
    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point, add adds them to neighAtt and dirDist. If forceFullCheck is false (the default) returns 
    /// as soon as the position is detected to be too far (without calculating the other distances).
    /// returns a DirDistances that contains the general distances (but no direction dependent data)
    DirDistances!(T) addDirsOf(DynPVector!(T,XType)newPos,Point newPoint,T pSize,ref LocalGrowableArray!(PointAndDir) neighAtt,
        ref LocalGrowableArray!(DirDistances!(T)) dirDist,bool forceFullCheck=false)
    {
        void log(void delegate(CharSink s)writer){
            logMsg(delegate void(CharSink sink){
                dumper(sink)(".addDirsOf")(newPoint.data)(" ");
                writer(sink);
            });
        }
        DirDistances!(T) dDist;
        if (newPoint==point && !forceFullCheck) {
            dDist.xDist=0;
            dDist.dualDist=0;
            dDist.cartesianDist=0;
            return dDist;
        }
        log(delegate void(CharSink s){s("starting");});
        bool veryClose=false;
        bool discarded=false;
        bool added=false;
        
        auto diff=newPos.dup();
        diff.opBypax(pos.dynVars.x,-1,1);
        
        auto internalDiff=diff.norm2(); // diff norm in internal coordinates
        dDist.xDist=internalDiff;
        if (dDist.veryFar(localContext)) {
            log((CharSink s){ dumper(s)("veryFar in internal coord"); });
            if (!forceFullCheck) return dDist;
        }
        auto deriv1=pos.dynVars.dVarStruct.emptyDx();
        deriv1[]=0;
        pos.addToTSpace!(T)(diff,deriv1);
        diff.giveBack();
        // make an optimized version when the overlap is the identity??
        auto deriv1Dual=pos.dynVars.dVarStruct.emptyDualDx();
        pos.toDualTSpace!(T)(deriv1,deriv1Dual);
        // original norm in dual space
        auto origNorm=deriv1Dual.norm2();
        dDist.dualDist=origNorm;
        log((CharSink s){ dumper(s)("distance is ")(dDist); });
        if (dDist.veryFar(localContext)) {
            log((CharSink s){ dumper(s)("veryFar"); });
            if (!forceFullCheck) return dDist;
        }
        if (origNorm<localContext.minNormDual()*explorationSize){
            log((CharSink s){ dumper(s)("veryClose"); });
            veryClose=true;
        }
        // norm in cartesian units (T space approx)
        auto normCartesian=sqrt(pos.dotInTSpace(deriv1,deriv1Dual));
        dDist.cartesianDist=normCartesian;
        log((CharSink s){ dumper(s)("did cartesian dist, now ")(dDist); });
        if (dDist.veryFar(localContext)) {
            log((CharSink s){ dumper(s)("veryFar"); });
            if (!forceFullCheck) return dDist;
        }
        if (normCartesian<localContext.minRealNormSelf2){
            log((CharSink s){ dumper(s)("veryClose"); });
            veryClose=true;
        }
        if (!dDist.neighbor(localContext,explorationSize,pSize)) {
            log((CharSink s){ dumper(s)("not neighbor"); });
            if (!forceFullCheck) return dDist;
        }
        bool fRef;
        synchronized(this) {
            fRef=hasFrameOfRef;
        }
        uint dirMax=invalidDir;
        if (fRef){
            log((CharSink s){ dumper(s)("checking frameOfRef")(" origNorm:")(origNorm)(" localContext.minNormDual")(localContext.minNormDual)(" explorationSize")(explorationSize)(" localContext.maxNormDual")(localContext.maxNormDual); });
            size_t iMax;
            bool didCalcDualDir=false;
            
            // check for ill defined directions in cartesian space
            auto deriv2=rotateVEi(minDir.toGroup!(DxType)(),0,deriv1);
            // length of the directions in cartesian units
            auto ov=pos.maybeDerivOverlap();
            NArray!(T,1) lenDirs;
            if (ov!is null){
                DynPMatrix!(T,DxType,DualDxType) rotOv;
                rotOv=new DynPMatrix!(T,DxType,DualDxType)(ov.colGroup,ov.rowGroup,ov.data.T.dup(true));
                //rotateVEi(minDir.toGroup!(DxType)(),0,rotOv);
                // should probably transpose and rotate again??? to check
                assert(0,"to do: should be checked");
                lenDirs=diag(rotOv.data);
            } else {
                lenDirs=ones!(T)(ndim); // should maybe avoid allocating so much using repeat...
            }
            unaryOpStr!("*aPtr0=sqrt(*aPtr0);")(lenDirs);
            
            void calcDualDir(){ 
                didCalcDualDir=true;
                // this direction is the "main" direction and should be the first
                // it is a neighbor in some direction
                DynPVector!(T,DualDxType) deriv2Dual;
                /// checking which direction of the dual space we are really exploring
                assert(!minDir.isDummy);
                deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
                dDist.minDirProj=deriv2Dual[0];
                volatile T maxVal=0;
                iMax=0;
                foreach (i,v;deriv2Dual.sLoop){
                    auto vp=abs(v);
                    if (vp>maxVal){
                        synchronized{
                            if (vp>maxVal){
                                maxVal=vp;
                                iMax=i;
                            }
                        }
                    }
                }
                dirMax=toDir(iMax,deriv2Dual[iMax]<0);
                log((CharSink s){ dumper(s)("mostly in the direction ")(dirMax)(", dimension ")(iMax); });
                bool noDir=false;
                auto cosVal=maxVal/origNorm;
                if (cosVal>localContext.sameDirCosAngle()){
                    log((CharSink s){ dumper(s)("mainly in same dir"); });
                    auto dOpt2=pow2(explorationSize-maxVal)+pow2(origNorm)-pow2(maxVal);
                    if (dOpt2<localContext.dirDualSize2){
                        log((CharSink s){ dumper(s)("within perfect dir radius"); });
                        auto dOpt=((dOpt2>0)?cast(T)sqrt(dOpt2):cast(T)0);
                        DirDistances!(T) dDist1=dDist;
                        neighAtt(PointAndDir(newPoint,dirMax));
                        added=true;
                        dDist1.dualDirDist=dOpt;
                        dDist1.dualCos=cosVal;
                        
                        // complete data of dDist
                        auto idealLen=explorationSize*lenDirs[iMax];
                        dDist1.cartesianCos=deriv2[iMax]/normCartesian;
                        auto dOpt2C=pow2(normCartesian)-pow2(deriv2[iMax])
                            +localContext.inDirCartesianScale2()*pow2(abs(deriv2[iMax])-idealLen);
                        dDist1.cartesianDirDist=dOpt2C;
                        
                        dirDist(dDist1);
                        auto actualDirFlags=exploredDirs.atomicCAS(dirMax,DirFlags.Explored,DirFlags.Free);
                        log((CharSink s){ dumper(s)("added direction"); });
                    }
                }
                deriv2Dual.giveBack();
            }
            
            if ((origNorm>localContext.minNormDual*explorationSize && origNorm<localContext.maxNormDual*explorationSize)||forceFullCheck){
                calcDualDir();
            }
            // loop on all directions checking cartesian thresholds
            auto dim=ndim;
            for(size_t idim=0;idim<dim;++idim){
                auto idealLen=explorationSize*lenDirs[idim];
                if (abs(normCartesian-idealLen)<localContext.minRealNormSelf2 &&
                    abs(deriv2[idim])>=localContext.sameDirCosAngleCartesian()*explorationSize*lenDirs[idim])
                {
                    auto dOpt2=pow2(normCartesian)-pow2(deriv2[idim])
                        +localContext.inDirCartesianScale2()*pow2(abs(deriv2[idim])-idealLen);
                    if (dOpt2<=localContext.dirCartesianSize2){
                        auto dirAtt=toDir(idim,deriv2[idim]<0);
                        if (dirAtt!=dirMax){
                            log((CharSink s){ dumper(s)("found almost linear dependent direction ")(dirAtt); });
                            exploredDirs.atomicCAS(dirAtt,DirFlags.Explored,DirFlags.Free);
                            if (idealLen<=localContext.zeroLen()){ // don't store directions for invalid dirs
                                log((CharSink s){ dumper(s)("avoiding storing linear dependent direction ")(dirAtt); });
                            } else {
                                if (!didCalcDualDir) calcDualDir();
                                DirDistances!(T) dDist1=dDist;
                                dDist1.cartesianCos=deriv2[idim]/normCartesian;
                                dDist1.cartesianDirDist=((dOpt2>0)?sqrt(dOpt2):0);
                                added=true;
                                neighAtt(PointAndDir(newPoint,dirAtt));
                                dirDist(dDist1);
                                log((CharSink s){ dumper(s)("storing almost linear dependent direction ")(dirAtt); });
                            }
                        }
                    }
                }
            }
        }
        if (veryClose || origNorm<=localContext.minNormDual*explorationSize){
            neighAtt(PointAndDir(newPoint,0));
            dirDist(dDist);
            if ((!added) && dirMax!=invalidDir) {
                neighAtt(PointAndDir(newPoint,dirMax));
                dirDist(dDist);
            }
            added=true;
            log((CharSink s){ dumper(s)("adding very close direction"); });
        }
        if ((!added)){
            if (origNorm<=localContext.maxNormDual*explorationSize+localContext.zeroLen()){
                added=true;
                neighAtt(PointAndDir(newPoint,ndirs)); // not in a specific direction
                dirDist(dDist);
                if (dirMax!=invalidDir) {
                    neighAtt(PointAndDir(newPoint,dirMax));
                    dirDist(dDist);
                }
                log((CharSink s){ dumper(s)("storing unspecific direction"); });
            } else if (origNorm<=localContext.maxNormDual*pSize+localContext.zeroLen()) {
                added=true;
                neighAtt(PointAndDir(newPoint,invalidDir)); // a neighbor of the other point
                dirDist(dDist);
                if (dirMax!=invalidDir) {
                    neighAtt(PointAndDir(newPoint,dirMax));
                    dirDist(dDist);
                }
                log((CharSink s){ dumper(s)("storing as neighbor of the other point only"); });
            } else {
                log((CharSink s){ dumper(s)("not a neighbor"); });
            }
        }
        log((CharSink s){ dumper(s)("final dDist:")(dDist); });
        return dDist;
    }
    /// adds the given point as neighbor
    void addNeighbor(Point p){
        synchronized(this){
            if (findFirstPred(neighbors.data,delegate bool(PointAndDir pDir){ return pDir.point==p; })!=neighbors.length){
                logMsg(delegate void(CharSink s){
                    dumper(s)("addNeighbor(")(p.data)(") already added\n");
                });
                return;
            }
        }
        auto tNow=ev_time();
        auto mPoint=localContext.createLocalPoint(p,tNow);
        scope(exit){ mPoint.release(); }
        DynPVector!(T,XType) newPos=mPoint.pos.dynVars.x;
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);
        DirDistances!(T)[128] bufDirDist;
        auto dirDist=lGrowableArray(bufDirDist,0);
        scope(exit){
            neighAtt.deallocData();
            dirDist.deallocData();
        }
        bool hasGrad=hasFrameOfRef;
        auto res=addDirsOf(newPos,p,mPoint.explorationSize,neighAtt,dirDist);
        if (neighAtt.length==0){ // should not be needed
            neighAtt(PointAndDir(p,invalidDir));
            dirDist(res);
        }
        logMsg(delegate void(CharSink s){
            dumper(s)("addNeighbor(")(p)(") will add ")(dirDist.data);
        });
        addNeighbors(neighAtt.data,dirDist.data,hasGrad);
    }
    /// evaluates with the given context, returns if an evaluate was really done
    bool evalWithContext(LocalCalculationContext c,bool alwaysGrad=false){
        version(TrackPNet){
            logMsg(delegate void(CharSink s){
                dumper(s)("evalWithContext(")(c.contextId)(",")(alwaysGrad)(")");
            });
        }
        bool calcE=(gFlags&GFlags.EnergyEvaluated)==0;
        bool calcF=((gFlags&GFlags.GradientEvaluated)==0 &&
         (localContext.gradEagerness()>=GradEagerness.Always||(gFlags&GFlags.GradientInProgress)!=0||alwaysGrad));
        if (calcE||calcF){
            Real e=pos.dynVars.potentialEnergy;
            Real eErr=pos.dynVars.potentialEnergyError;
            Real diff2=0;
            mixin(withPSys(`
            pSys.checkX();
            foreach(i,ref v;pSys.dynVars.x.sLoop){ /+ should do an optimized native impl +/
                auto x=pos.dynVars.x[i];
                diff2+=pow2(v-x);
                v=x;
            }
            //pSys.dynVars.x[]=pos.dynVars.x;
            pSys.dynVars.potentialEnergy=pos.dynVars.potentialEnergy;
            pSys.dynVars.potentialEnergyError=pos.dynVars.potentialEnergyError;
            `,"c."));
            c.changedDynVars(ChangeLevel.SmallPosChange,sqrt(diff2));
            char[64] buf;
            auto arr=lGrowableArray(buf,0,GASharing.Local);
            writeOut(&arr.appendArr,point.data);
            auto extRef=arr.takeData;
            c.externalRefSet(extRef);
            version(TrackPNet) logMsg(delegate void(CharSink s){ dumper(s)("will calculate EF(")(calcE)(",")(calcF)(")"); });
            c.updateEF(calcE,calcF);
            version(TrackPNet) logMsg(delegate void(CharSink s){ s("did calculate EF"); });
            if (calcE) {
                pos.dynVars.potentialEnergy=c.potentialEnergy();
                pos.dynVars.potentialEnergyError=c.potentialEnergyError();
                if ((!isNaN(e))&&(e!=pos.dynVars.potentialEnergy)){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("energy of point ")(point.data)(" changed from ")(e)(" to ")(pos.dynVars.potentialEnergy);
                    }),__FILE__,__LINE__);
                }
            }
            if (calcF) {
                pos.checkMddx();
                mixin(withPSys(`
                pos.dynVars.mddx[]=pSys.dynVars.mddx;
                pos.dynVars.mddxError=pSys.dynVars.mddxError;`,"c."));
            }
            e=pos.dynVars.potentialEnergy;
            eErr=pos.dynVars.potentialEnergyError;
            if (isNaN(e)){
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy after calculation in point ")(point.data);
                }),__FILE__,__LINE__);
            }
            auto target=localContext.ownerOfPoint(point);
            bool calcMore=false;
            if (calcE && (!calcF)){
                if (localContext.gradEagerness()==GradEagerness.Speculative){
                    version(TrackPNet) logMsg(delegate void(CharSink s){ s("ask for extra gradient"); });
                    calcMore=localContext.speculativeGradientLocal(target,point,e,eErr);
                } else {
                    version(TrackPNet) logMsg(delegate void(CharSink s){ s("will communicate energy"); });
                    localContext.addEnergyEvalLocal(target,point,e,eErr);
                }
            }
            if (calcMore){
                version(TrackPNet) logMsg(delegate void(CharSink s){ s("calc extra gradient"); });
                calcF=true;
                pos.checkMddx();
                c.updateEF(false,true);
                mixin(withPSys(`
                pos.dynVars.mddx[]=pSys.dynVars.mddx;
                pos.dynVars.mddxError=pSys.dynVars.mddxError;`,"c."));
                // update also the energy? it should not change...
            }
            if (calcF){
                version(TrackPNet) logMsg(delegate void(CharSink s){ s("communicate force eval"); });
                if (isLocalCopy){
                    localContext.addGradEvalLocal(target,point,pSysWriter(pos));
                } else {
                    didGradEval();
                }
            }
            version(TrackPNet){
                logMsg(delegate void(CharSink s){
                    dumper(s)(" did evalWithContext(")(c.contextId)(",")(alwaysGrad)(")");
                });
            }
            return true;
        } else {
            version(TrackPNet){
                logMsg(delegate void(CharSink s){
                    dumper(s)(" no calculation in evalWithContext(")(c.contextId)(",")(alwaysGrad)(")");
                });
            }
            return false;
        }
    }
    /// notifies a GFlags change from the given gFlags
    /// sends the address of a stack allocated GFlagsChange, avoids allocation, but you have to take care
    /// as it won't remain valid after the notification
    void notifyGFlagChange(uint oldGFlags){
        if (gFlags!=oldGFlags){
            GFlagsChange gFlagsChange;
            gFlagsChange.oldGFlags=oldGFlags;
            gFlagsChange.newGFlags=gFlags;
            gFlagsChange.point=point;
            if ((gFlags&GFlags.LocalCopy)==0){
                localContext.nCenter.notify("localPointChangedGFlags",box(&gFlagsChange));
            } else {
                localContext.nCenter.notify("copiedPointChangedGFlags",box(&gFlagsChange));
            }
        }
    }
    /// notifies that the attractor did change (could be optimized a bit to ensure that double notifications are avoided (keeping id of bcasted))
    void notifyAttractorChange(){
        PointAndDir[128] buf;
        Real[128] buf2;
        auto lMem=LocalMem(buf2);
        auto localNeigh=lGrowableArray(buf,0);
        PointEMin eAndMin;
        synchronized(this){
            localNeigh.appendArr(neighbors.data);
            eAndMin.point=point;
            eAndMin.minimum=_attractor.minimum;
            eAndMin.id=_attractor.id;
            eAndMin.energy=energy;
        }
        auto lData=localNeigh.data;
        auto nPts=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return lData[i].point; },
            lData.length);
        localNeigh.deallocData();
        foreach (i;nPts.pLoop){
            auto silo=i.owner;
            localContext.neighborHasEnergy(silo,i.points,eAndMin);
        }
        localContext.nCenter.notify("attractorChange",box(&eAndMin));
    }
    /// the energy of a neighbor was calculated (or the minimum updated)
    void localNeighEnergy(PointAndDir[] neigh,DirDistances!(T)[] dirDist,PointEMin eAndMin){
        typeof(_gFlags) newFlags;
        synchronized(this){
            if (isNaN(energy) || (gFlags&GFlags.EnergyEvaluated)==0) {
                /+assert((gFlags&GFlags.EnergyInfo)==GFlags.InProgress,collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy and not in progress in point ")(this); }));+/
                logMsg(delegate void(CharSink s){
                    dumper(s)("skipping energy info (")(eAndMin)(") of neighbor ")(neigh)
                        (" as this point doesn't have a calculated energy yet");
                });
                return;
            }
        }
        if (dirDist.length>0 && dirDist[0].energy==energy&&_attractor.throughPoint!=eAndMin.point) return; // quick return for non relevant minimum updates
        auto e=eAndMin.energy;
        if (e<energy){
            newFlags|=GFlags.NeighValsDecrease;
        }else if (e>energy){
            newFlags|=GFlags.NeighValsIncrease;
        } else {
            newFlags|=GFlags.NeighValsSame;
        }
        // check for special directions
        uint lastDir=ndirs-1;
        int specialDir=0;
        bool forceDir=false;
        for(size_t iNeigh=0;iNeigh<neigh.length;++iNeigh){
            if (neigh[iNeigh].dir==1){ // along the force
                specialDir=1;
                forceDir=true;
                if (e<=energy) newFlags|=GFlags.AlongForcesDecrease;
                if (e>=energy) newFlags|=GFlags.AlongForcesIncrease;
            }
            if (neigh[iNeigh].dir==lastDir){ // along the gradient
                specialDir=1;
                if (e<=energy) newFlags|=GFlags.AlongGradientDecrease;
                if (e>=energy) newFlags|=GFlags.AlongGradientIncrease;
            } else if (neigh[iNeigh].dir>1 && neigh[iNeigh].dir<lastDir && specialDir==0){
                specialDir=-1;
            }
            if (dirDist.length>iNeigh) dirDist[iNeigh].energy=e;
        }
        if (specialDir==-1&&dirDist.length>0 && (!isNaN(dirDist[0].minDirProj))){
            auto eCorr=e;
            eCorr+=dirDist[0].minDirProj*minDirScale();
            if (eCorr<energy && (gFlags&GFlags.CorrNeighValsDecrease)==0){
                newFlags|=GFlags.CorrNeighValsDecrease;
            }else if (eCorr>energy){
                newFlags|=GFlags.CorrNeighValsIncrease;
            } else {
                newFlags|=GFlags.CorrNeighValsSame;
            }
        }
        if ((forceDir||specialDir==-1||(!hasFrameOfRef)) && (isNaN(_attractor.energyThroughPoint) || _attractor.energyThroughPoint>e||
            (_attractor.throughPoint==eAndMin.point && _attractor.idThroughPoint<eAndMin.id)))
        {
            bool didChangeMinimum=false;
            synchronized(this){
                if ((forceDir||specialDir==-1||(!hasFrameOfRef))&&(isNaN(_attractor.energyThroughPoint) || _attractor.energyThroughPoint>e || 
                    (_attractor.throughPoint==eAndMin.point && _attractor.idThroughPoint<eAndMin.id)))
                {
                    _attractor.throughPoint=eAndMin.point;
                    _attractor.energyThroughPoint=eAndMin.energy;
                    _attractor.idThroughPoint=eAndMin.id;
                    if (_attractor.minimum!=eAndMin.minimum){
                        didChangeMinimum=true;
                        _attractor.minimum=eAndMin.minimum;
                        ++(_attractor.id);
                    }
                }
            }
            if (didChangeMinimum) {
                notifyAttractorChange();
            }
        }
        if ((gFlags|newFlags)!=gFlags){
            uint oldGFlags,newGFlags;
            synchronized(this){
                oldGFlags=gFlags;
                _gFlags|=newFlags;
                newGFlags=_gFlags;
            }
            if (oldGFlags!=newGFlags) notifyGFlagChange(oldGFlags);
        }
    }
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    /// it is possible that only the minimum has changed
    void addEnergyEvalOther(PointEMin energyOther){
        auto p=energyOther.point;
        PointAndDir[128] buf;
        DirDistances!(T)[64] buf2;
        LocalGrowableArray!(PointAndDir) localCopy=lGrowableArray(buf,0);
        LocalGrowableArray!(DirDistances!(T)) localCopyDir=lGrowableArray(buf2,0);
        // find neighbor
        synchronized(this){
            size_t iStart=0;
            while (1){
                while(iStart<neighbors.length){
                    if (p==_neighbors[iStart].point) break;
                    ++iStart;
                }
                if (iStart==neighbors.length) break;
                size_t iEnd=iStart+1;
                while(iEnd<neighbors.length){
                    if (p!=_neighbors[iEnd].point) break;
                    ++iEnd;
                }
                localCopy.appendArr(_neighbors.data[iStart..iEnd]);
                if (_neighDistances.length>=iEnd){
                    localCopyDir.appendArr(_neighDistances.data[iStart..iEnd]);
                }
                iStart=iEnd;
            }
        }
        if (localCopy.length!=0){
            localNeighEnergy(localCopy.data,localCopyDir.data,energyOther);
        }
        localCopy.deallocData();
        localCopyDir.deallocData();
    }
    
    /// if the gradient should be calculated, if returns true, it assumes that the gradient is then in progress
    bool shouldCalculateGradient(){
        assert(!isLocalCopy);
        bool res=false;
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags;
            if ((gFlags&GFlags.GradientInfo)==0){
                _gFlags=(gFlags & ~GFlags.GradientInfo)|GFlags.GradientInProgress;
                exploredDirs.atomicCAS(0,DirFlags.Explored,DirFlags.Free);
                return true;
            }
        }
        notifyGFlagChange(oldGFlags);
        return false;
    }
    
    /// sets the energy
    /// is really perfomed only the first time (and returns true)
    /// if true is returned you should call didEnergyEval
    bool setEnergy(Real e,Real eError){
        assert(!isLocalCopy);
        if (isNaN(e)){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("energy should not be NAN in point ")(point);
            }),__FILE__,__LINE__);
        }
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags;
            switch(oldGFlags&GFlags.EnergyInfo){
            case GFlags.None, GFlags.EnergyInProgress:
                pos.dynVars.potentialEnergy=e;
                pos.dynVars.potentialEnergyError=eError;
                oldGFlags=gFlags;
                _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
                break;
            case GFlags.EnergyKnown,GFlags.EnergyEvaluated:
                if (e!=pos.dynVars.potentialEnergy){
                    logMsg(delegate void(CharSink s){
                        dumper(s)("Warning: double energy eval with different value ignored:")
                            (e)(" vs ")(pos.dynVars.potentialEnergy);
                    });
                }
                return false;
            default:
                assert(0);
            }
        }
        notifyGFlagChange(oldGFlags);
        return true;
    }
    /// energy must be valid, should be called just once
    /// gFlags is updated at the end. Sends (and gets) the energy to all neighbors
    void didEnergyEval(){
        assert(!isLocalCopy);
        switch(gFlags&GFlags.EnergyInfo){
        case GFlags.None, GFlags.EnergyInProgress:
            throw new Exception("energy should be known when calling didEnergyEval",__FILE__,__LINE__);
        case GFlags.EnergyKnown:
            break;
        case GFlags.EnergyEvaluated:
            throw new Exception("energy should not be already broadcasted when calling didEnergyEval",__FILE__,__LINE__);
        default:
            assert(0);
        }
        updateNeighE();
        auto s=owner();
        localContext.notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
            obs.addEnergyEvalLocal(s,point,energy,energyError);
        });
    }
    /// gets the energies of the neighbors and broadcasts its own energy to them
    void updateNeighE(){
        PointAndDir[128] buf;
        Real[128] buf2;
        auto lMem=LocalMem(buf2);
        auto localNeigh=lGrowableArray(buf,0);
        synchronized(this){
            localNeigh.appendArr(neighbors.data);
        }
        auto lData=localNeigh.data;
        auto nPts=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return lData[i].point; },
            lData.length);
        localNeigh.clearData();
        auto energies=lMem.allocArr!(PointEMin)(nPts.points.length);
        uint newGFlags=0;
        foreach (i;nPts.pLoop){
            auto silo=i.owner;
            auto energ=localContext.energyForPointsLocal(silo,i.points,energies[i.localLb..i.localUb]);
            auto nPts=i.points.length;
            for(size_t iPts=0;iPts<nPts;++iPts){
                addEnergyEvalOther(energ[i.localLb+iPts]);
            }
        }
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags;
            _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyEvaluated;
            if ((_gFlags&GFlags.NeighValsDecrease)==0){
                assert(isNaN(_attractor.energyThroughPoint)||_attractor.energyThroughPoint>=energy);
                _attractor.energyThroughPoint=energy;
                _attractor.minimum=point;
                _attractor.throughPoint=point;
                ++(_attractor.id);
                _attractor.idThroughPoint=_attractor.id;
            }
        }
        notifyGFlagChange(oldGFlags);
        notifyAttractorChange();
    }
    
    /// the gradient of a neighbor was calculated (and possibly also the energy for the first time)
    void localNeighGrad(PointAndDir[] neigh,DirDistances!(T)[] dirDist,LazyMPLoader!(T)mainPoint,PointEMin eAndMin){
        typeof(_gFlags) newFlags;
        synchronized(this){
            if (isNaN(energy) || (gFlags&GFlags.EnergyEvaluated)==0) {
                /+assert((gFlags&GFlags.EnergyInfo)==GFlags.InProgress,collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy and not in progress in point ")(this); }));+/
                logMsg(delegate void(CharSink s){
                    dumper(s)("skipping energy info (")(eAndMin)(") of neighbor ")(neigh)
                        (" as this point doesn't have a calculated energy yet");
                });
                return;
            }
        }
        auto e=eAndMin.energy;
        if (e<energy){
            newFlags|=GFlags.NeighValsDecrease;
        } else if (e>energy){
            newFlags|=GFlags.NeighValsIncrease;
        } else {
            newFlags|=GFlags.NeighValsSame;
        }
        // check for special directions
        uint lastDir=ndirs-1;
        int specialDir=0;
        bool forceDir=false;
        for(size_t iNeigh=0;iNeigh<neigh.length;++iNeigh){
            if (neigh[iNeigh].dir==1){ // along the force
                specialDir=1;
                forceDir=true;
                if (e<=energy) newFlags|=GFlags.AlongForcesDecrease;
                if (e>=energy) newFlags|=GFlags.AlongForcesIncrease;
                if (hasFrameOfRef){
                    // to do: this is the place to put gradient based topology info, will be useful with free energy
                    //auto nPoint=mainPoint.mainPoint(localContext);
                    //nPoint.pos.dynVars.mddx // scalar prod with minDir, expected to be >0 (if it is a "normal" point), comparing with minDirScale gives an idea of the curvature
                }
            }
            if (neigh[iNeigh].dir==lastDir){ // along the gradient
                specialDir=1;
                if (e<=energy) newFlags|=GFlags.AlongGradientDecrease;
                if (e>=energy) newFlags|=GFlags.AlongGradientIncrease;
                if (hasFrameOfRef){
                    // to do: this is the place to put gradient based topology info, will be useful with free energy
                    //auto nPoint=mainPoint.mainPoint(localContext);
                    //nPoint.pos.dynVars.mddx // scalar prod with minDir, expected to be >0 (if it is a "normal" point), comparing with minDirScale gives an idea of the curvature
                }
            } else if (neigh[iNeigh].dir>1 && neigh[iNeigh].dir<lastDir && specialDir==0){
                specialDir=-1;
            }
            if (dirDist.length>iNeigh) dirDist[iNeigh].energy=e;
        }
        if (specialDir==-1&&dirDist.length>0 && (!isNaN(dirDist[0].minDirProj))){
            auto eCorr=e;
            eCorr+=dirDist[0].minDirProj*minDirScale();
            if (eCorr<energy && (gFlags&GFlags.CorrNeighValsDecrease)==0){
                newFlags|=GFlags.CorrNeighValsDecrease;
            }else if (eCorr>energy){
                newFlags|=GFlags.CorrNeighValsIncrease;
            } else {
                newFlags|=GFlags.CorrNeighValsSame;
            }
        }
        if ((forceDir||specialDir==-1) && (isNaN(_attractor.energyThroughPoint) || _attractor.energyThroughPoint>e||
            (_attractor.throughPoint==eAndMin.point && _attractor.idThroughPoint<eAndMin.id)))
        {
            bool changedAttractor=false;
            synchronized(this){
                if (isNaN(_attractor.energyThroughPoint) || _attractor.energyThroughPoint>e||
                    (_attractor.throughPoint==eAndMin.point && _attractor.idThroughPoint<eAndMin.id))
                {
                    _attractor.throughPoint=eAndMin.point;
                    _attractor.energyThroughPoint=eAndMin.energy;
                    _attractor.idThroughPoint=eAndMin.id;
                    if (_attractor.minimum!=eAndMin.minimum){
                        _attractor.minimum=eAndMin.minimum;
                        ++(_attractor.id);
                        changedAttractor=true;
                    }
                }
            }
            if (changedAttractor){
                notifyAttractorChange();
            }
        }
        if ((gFlags|newFlags)!=gFlags){
            uint oldGFlags,newGFlags;
            synchronized(this){
                oldGFlags=gFlags;
                _gFlags|=newFlags;
                newGFlags=_gFlags;
            }
            if (oldGFlags!=newGFlags) notifyGFlagChange(oldGFlags);
        }
    }
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    void addGradEvalOther(LazyMPLoader!(T)mp,PointEMin eAndMin){
        PointAndDir[128] buf;
        DirDistances!(T)[64] buf2;
        LocalGrowableArray!(PointAndDir) localCopy=lGrowableArray(buf,0);
        LocalGrowableArray!(DirDistances!(T)) localCopyDir=lGrowableArray(buf2,0);
        Point p=mp.point;
        // find neighbor
        synchronized(this){
            size_t iStart=0;
            while (1){
                while(iStart<neighbors.length){
                    if (p==_neighbors[iStart].point) break;
                    ++iStart;
                }
                if (iStart==neighbors.length) break;
                size_t iEnd=iStart+1;
                while(iEnd<neighbors.length){
                    if (p!=(*(neighbors()))[iEnd].point) break;
                    ++iEnd;
                }
                localCopy.appendArr(neighbors.data[iStart..iEnd]);
                if (_neighDistances.length>=iEnd)
                    localCopyDir.appendArr(_neighDistances.data[iStart..iEnd]);
                iStart=iEnd;
            }
        }
        if (localCopy.length!=0){
            localNeighGrad(localCopy.data,localCopyDir.data,mp,eAndMin);
        }
        localCopy.deallocData();
    }
    /// marks this point as dropped returns true if it hadn't been dropped yet...
    bool drop(){
        synchronized(this){
            auto oldFlags=_gFlags;
            if (oldFlags&GFlags.LocalCopy){
                throw new Exception("cannot drop local copies of point",__FILE__,__LINE__);
            }
            _gFlags|=GFlags.OldApprox|GFlags.DoNotExplore;
            if ((oldFlags&GFlags.OldApprox)==0){
                return true;
            }
            return false;
        }
    }
    /// adds the gradient to this point
    /// energy must be valid
    /// gFlags is updated at the end assigning directions to all neighbors
    /// attractor might be modified
    void didGradEval(){
        assert(!isLocalCopy);
        version (TrackPNet) logMsg(delegate void(CharSink s){ s("didGradEval"); });
        bool newE=false;
        uint oldGFlags,newGFlags;
        if ((gFlags&GFlags.HasRefFrame)==0){
            buildMinDir();
        }
        PointAndDir[128] buf;
        auto localNeigh=lGrowableArray(buf,0);
        DirDistances!(T)[64] buf2;
        auto localNeighDir=lGrowableArray(buf2,0);
        synchronized(this){
            oldGFlags=gFlags;
            _gFlags|=GFlags.HasRefFrame;
            if ((gFlags&GFlags.EnergyEvaluated)==0){
                _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
                if (isNaN(energy())) {
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("invalid energy after evaluation in point ")(point);
                    }),__FILE__,__LINE__);
                }
                newE=true;
            }
            if ((oldGFlags & GFlags.GradientEvaluated)!=0) {
                logMsg(delegate void(CharSink s){ 
                    s("Warning, double eval of gradient, if the gradient changed this means inconsistent minDir and gradient.");
                    if (!newE) s("\nExiting didGradEval now.");
                });
                if (!newE) return;
            } else {
                exploredDirs.atomicCAS(0,DirFlags.Explored,DirFlags.Free);
                _gFlags=(gFlags & ~GFlags.GradientInfo)|GFlags.GradientKnown;
            }
            newGFlags=_gFlags;
            localNeigh.appendArr(neighbors.data);
            localNeighDir.appendArr(neighDistances.data);
            neighbors.clearData();
            neighDistances.clearData();
            // clear attractor (would be better to clear it only if it is in the "wrong" direction, change??)
            if (_attractor.minimum !=  point){
                _attractor.energyThroughPoint=Real.init;
                _attractor.minimum=Point(0);
                _attractor.throughPoint=Point(0);
                ++(_attractor.id);
                _attractor.idThroughPoint=0;
            }
        }
        if (newGFlags!=oldGFlags){
            notifyGFlagChange(oldGFlags);
        }
        {
            auto lData=localNeigh.data;
            auto nPts=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return lData[i].point; },
                lData.length);
            Real[] energies;
            auto tNow=ev_time();
            foreach (i;nPts.pLoop){
                auto nPts=i.points.length;
                for(size_t iPts=0;iPts<nPts;++iPts){
                    auto pointAtt=i.points()[iPts];
                    auto pAtt=localContext.createLocalPoint(pointAtt,tNow);
                    scope(exit) { pAtt.release(); }
                    checkIfNeighbor(pAtt.pos.dynVars.x, pointAtt,pAtt.explorationSize);
                    if (!isNaN(pAtt.energy)){
                        if (pAtt.hasFrameOfRef){
                            auto otherP=LazyMPLoader!(T)(pointAtt,true,tNow);
                            addGradEvalOther(otherP,pAtt.pointEMin);
                            otherP.release();
                        } else {
                            addEnergyEvalOther(pAtt.pointEMin);
                        }
                    }
                }
            }
        }
        auto ownr=owner;
        localContext.notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
            obs.addGradEvalLocal(ownr,point,pSysWriter(pos));
            if(newE) // we do it after to help when playing back the journal (skip???)
                obs.addEnergyEvalLocal(ownr,point,energy,energyError);
        });
        localNeigh.clearData();
        localNeighDir.clearData();
        PointEMin eAndMin;
        synchronized(this){
            oldGFlags=_gFlags;
            _gFlags=(gFlags & ~(GFlags.EnergyInfo|GFlags.GradientInfo))|
                (GFlags.EnergyEvaluated|GFlags.GradientEvaluated);
            newGFlags=_gFlags;
            if ((_gFlags&GFlags.NeighValsDecrease)==0){
                if (_attractor.minimum!=point){
                    assert(isNaN(_attractor.energyThroughPoint)||_attractor.energyThroughPoint<=energy);
                    _attractor.minimum=point;
                    _attractor.throughPoint=point;
                    _attractor.energyThroughPoint=energy;
                    ++(_attractor.id);
                    _attractor.idThroughPoint=_attractor.id;
                }
            }
            eAndMin.point=point;
            eAndMin.minimum=_attractor.minimum;
            eAndMin.id=_attractor.id;
            eAndMin.energy=energy;
            localNeigh.appendArr(neighbors.data);
            localNeighDir.appendArr(neighDistances.data);
        }
        if (oldGFlags!=newGFlags) {
            notifyGFlagChange(oldGFlags); // notify also earlier with GradientKnown???
        }
        {
            auto lData=localNeigh.data;
            auto tNow=ev_time();
            auto thisP=LazyMPLoader!(T)(point,true,tNow);
            scope nPts=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return lData[i].point; },
                lData.length);
            foreach (i;nPts.pLoop){
                auto silo=i.owner;
                localContext.neighborHasGradient(silo,thisP,i.points,eAndMin);
            }
            
        }
    }
    /// builds the direction toward the minimum
    void buildMinDir(){
        assert(!isNullT(pos.dynVars.mddx),"needs gradient");
        auto newMinDir=pos.dynVars.dVarStruct.emptyDualDx();
        auto mddxConstr=pos.dynVars.mddx.dup();
        localContext.constraints().applyDR(pos,mddxConstr);
        pos.toDualTSpace(mddxConstr,newMinDir);
        mddxConstr.giveBack();
        auto n2=newMinDir.norm22();
        if (n2==0){
            newMinDir[0]=1;
            _minDirScale=0;
        } else {
            _minDirScale=sqrt(n2);
            newMinDir*=cast(T)(1/minDirScale);
        }
        synchronized(this){
            _minDir=newMinDir;
        }
    }
    
    void retain(){
        if (atomicAdd(refCount,cast(size_t)1)==0){
            throw new Exception("refCount was 0 in retain",__FILE__,__LINE__);
        }
    }
    void release(){
        synchronized(localContext){
            auto oldV=atomicAdd(refCount,-cast(size_t)1);
            if (oldV==0){
                throw new Exception("refCount was 0 before release",__FILE__,__LINE__);
            }
            if (oldV==1){
                if ((gFlags&GFlags.LocalCopy)!=0){
                    localContext.dropCachedPoint(this);
                }
            }
        }
    }
    void opSliceAssign(DriedPoint!(T)p){
        _point=p.point;
        p.pos.copyTo(pos);
        if (!p.minDir.isDummy){
            if (_minDir.isDummy){
                _minDir=pos.dynVars.dVarStruct.emptyDualDx();
            }
            _minDir[]=p.minDir;
        }
        _exploredDirs=p.exploredDirs;
        _neighbors=p.neighbors;
        _attractor=p.attractor;
        _minDirScale=p.minDirScale;
        _explorationSize=p.explorationSize;
        _gFlags=p.gFlags;
    }
}

private MainPoint!(Real) dummyP; // to get some errors already compiling this
