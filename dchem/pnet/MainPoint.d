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
import blip.sync.Atomic;
import blip.math.Math:max,abs,pow2;
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

    /// internal help structure (to reorder stuff)
    struct PIdxOwn{
        Point point;
        size_t idx;
        SKey owner;
    }
    /// returns a PointToOwner that contains the points extracted with getPnt ordered by owner and point
    /// if compressPts is true (default) then duplicate points are removed
    /// if keepIdx is true (the default) the indexes of the points in the origina array are kept
    this(SKey delegate(Point) ownerMap,Point delegate(size_t i) getPnt, size_t nPts, bool compressPts=true,bool keepIdx=true){
        if (nPts==0) return;

        PIdxOwn[128] buf;
        auto pts=lGrowableArray(buf,1);
        scope(exit){ pts.deallocData(); }
        size_t i=0,j=1;
        pts[0].point=getPnt(0);
        pts[0].idx=0;
        pts[0].owner=ownerMap(getPnt(0));
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
        sort(points,delegate bool(PIdxOwn a,PIdxOwn b){ return a.owner<b.owner || (a.owner==b.owner && a.point<b.point); });
        if (compressPts){
            points=compress(points,delegate bool(PIdxOwn a,PIdxOwn b){ return a.point==b.point; });
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
    mixin(serializeSome("dchem.GFlagsChange","oldGFlags|newGFlags|point"));
    mixin printOut!();
}

/// structure to describe new neighbors
struct PointNeighbors(T){
    Point point;
    PointAndDir[] neighbors;
    DirDistances!(T)[] dirDist;
    Time time;
    mixin(serializeSome("dchem.PointNeighbors!("~T.stringof~")","point|neighbors|dirDist"));
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
    Point attractor(){ return _attractor; }
    Point _attractor;
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
    LocalGrowableArray!(PointAndDir) neighbors(){ return _neighbors; }
    LocalGrowableArray!(PointAndDir) _neighbors;
    /// neighbor distances, stored only if requested. 6 numbers for each neighbor:
    /// dualDist, cartDist, cosDual, rDistDual, cartCos, cartRDist
    LocalGrowableArray!(DirDistances!(T)) neighDistances(){ return _neighDistances; }
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
    size_t refCount;
    PoolI!(MainPoint) pool;
    
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
        return (gFlags&GFlags.PointBcasted)!=0;
    }
    /// if the frame of reference is set
    bool hasFrameOfRef(){
        return (gFlags&GFlags.GradientInfo)==GFlags.GradientEvaluated;
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
        res.neighbors=neighbors;
        res.minDirScale=minDirScale;
        res.explorationSize=explorationSize;
        return res;
    }
    /// if given the neighbor are taken into ownership, and freed by this object.
    this(LocalSilosI!(T) localContext, PoolI!(MainPoint) p,PointAndDir[] neighbors=[]){
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
    }
    /// constructor
    this(LocalSilosI!(T) localContext,Point point,ParticleSys!(T) pos,uint gFlags,
        T explorationSize,DynPVector!(T,DualDxType) minDir,T minDirScale=T.init,
        FlagsArray exploredDirs=null,PointAndDir[] neighbors=null,PoolI!(MainPoint)pool=null)
    {
        this._localContext=localContext;
        this._point=point;
        this._pos=pos;
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
        synchronized(this){
            if (_gFlags==oldVal){
                _gFlags=newVal;
            }
            return _gFlags;
        }
    }
    /// atomic op on the flags of this point
    uint gFlagsAtomicOp(uint delegate(uint) op){
        synchronized(this){
            auto oldV=_gFlags;
            _gFlags=op(oldV);
            return oldV;
        }
    }
    /// explores the next direction, immediately marks it as in exploration
    /// returns direction 0 if the gradient has to be evaluated, an invalid direction if the current point is in evaluation,
    /// and an invalid point only if all direction are explored/ the point should not be explored
    /// if lastIsLast is true then the last direction (gradient) is explored as last
    PointAndDir exploreNext(FlagsArray methodDirs=null,bool lastIsLast=true){
        auto exploreFlags=DirFlags.Explored;
        if ((gFlags&GFlags.DoNotExplore)!=0){
            return PointAndDir(Point(0),invalidDir);
        }
        if ((gFlags&GFlags.GradientEvaluated)==0){
            synchronized(this){
                if ((gFlags&(GFlags.InProgress|GFlags.GradientEvaluated))==0){
                    _gFlags|=GFlags.GradientInProgress;
                    return PointAndDir(point,0);
                }
            }
        }
        if ((gFlags&GFlags.InProgress)!=0){
            return PointAndDir(this.point,uint.max); // needs to wait evaluation
        }
        if ((gFlags&GFlags.FullyExplored)==0){
            synchronized(this){
                if ((gFlags&GFlags.FullyExplored)==0){
                    size_t nextDir=exploredDirs.length;
                    if(lastIsLast && methodDirs.atomicCAS(1,exploreFlags,0)==0 && 
                        exploredDirs.atomicCAS(1,exploreFlags,0)==0){ // explore down first
                        nextDir=0;
                    } else {
                        auto start=localContext.rand.uniformR(nextDir);
                        nextDir=exploredDirs.findFreeAndSet(start,exploreFlags,lastIsLast,methodDirs);
                    }
                    if (nextDir!=exploredDirs.length) return PointAndDir(this.point,nextDir);
                    _gFlags|=GFlags.FullyExplored;
                }
            }
        }
        return PointAndDir(Point(0),0);
    }
    /// returns the number of dimensions
    uint ndim(){
        auto g=pos.dynVars.dVarStruct.dualDxGroup;
        return g.posStruct.length+g.orientStruct.length+g.dofStruct.length;
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
        if(dir==ndirs-1) dir=0;
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
    /// calculates the position exploring from here in the given direction
    /// returns null if the position is *really* too close
    ParticleSys!(T) createPosInDir(uint dir){
        auto newDir=createDualDir(dir);
        auto cartesianNorm=pos.projectInDualTSpace(newDir);
        auto projNorm=newDir.norm2();
        if (cartesianNorm<=localContext.minRealNormSelf0){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)(" discarded because cartesian norm is too small. Projected norm was ")(projNorm)(".");
            });
            return null; // continue all the same or return pos???
        }
        assert(projNorm>0,"projNorm should not be 0 without cartesian norm being 0 too");
        auto scaleAtt=explorationSize/projNorm; // always normalizing in the non null directions
        /// immediately check expected cartesian norm change and possibly increase length
        /// this could be wrong for periodic directions, but is most likely better than not doing it
        if (cartesianNorm*scaleAtt<localContext.minRealNormSelf1){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)(" has small cartesian norm, increasing search length");
            });
            scaleAtt=localContext.minRealNormSelf1/cartesianNorm;
        }
        newDir*=scaleAtt;
        
        auto resPSys=pos.dup(PSDupLevel.DynProperties|PSDupLevel.DynPNullDx|PSDupLevel.DynPNullMddx|PSDupLevel.HiddenVars);
        pos.addFromDualTSpace!(T)(newDir,resPSys.dynVars.x);
        newDir.giveBack();
        auto err=localContext.constraints.applyR(resPSys);
        if (err>0.1){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)("discarded because we might have a constraint clash");
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
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in the internal coordinates of exploration direction ")(dir)(" of point ")
                    (point)(" is too small in internal coordinates ")(mainDists.xDist)(", discarding evaluation");
            });
            return false;
        }
        // original norm in dual space
        auto origNorm=mainDists.dualDist;
        if (origNorm<=localContext.minNormDualSelf()*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is too small:")(origNorm)(", discarding evaluation");
            });
            return false;
        }
        if (origNorm>localContext.maxNormDual*explorationSize || origNorm<localContext.minNormDual*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is outside the exploration zone with length ")(origNorm/explorationSize)("*")(explorationSize)
                    (", continuing evaluation");
            });
        }
        // norm in cartesian units
        auto normCartesian=mainDists.cartesianDist;
        if (normCartesian<localContext.minRealNormSelf2){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in real space of exploration in direction ")(dir)(" of point ")(point)
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
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("exploration in direction ")(dir)(" of point ")(point)
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
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is  mostly in the expected direction in the dual space, but ")
                                    (iMax)(", but not in the expected region, continuing evaluation");
                            });
                            nonVisited=true;
                        } else {
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is  very close to the starting point, and not in the expected direction");
                            });
                        }
                    } else {
                        auto actualDirFlags=exploredDirs.atomicCAS(iMax,DirFlags.Explored,DirFlags.Free);
                        if (actualDirFlags == DirFlags.Free){
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is not mostly in the expected direction in the dual space, but in dir ")
                                    (neighAtt[idir].dir)(", continuing evaluation declaring also it as visited");
                            });
                            nonVisited=true;
                        } else {
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is not mostly in the expected direction in the dual space, but in dir ")
                                    (neighAtt[idir].dir)(", and that was already visited");
                            });
                        }
                    }
                }
                if (!nonVisited) {
                    localContext.logMsg(delegate void(CharSink s){
                        dumper(s)("exploration in direction ")(dir)(" of point ")(point)
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
        if (neighs.length==0) return true; // return false?
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
            if (!hadGradient && hasFrameOfRef){
                gradChanged=true;
            } else {
                neighbors.appendArr(neighs); // duplicates are stored contiguously
            }
        }
        if (gradChanged){
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
        auto neigToNotify=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return neighs[i].point; },neighs.length);
        foreach (iPts;neigToNotify.pLoop){
            localContext.addPointToLocalNeighs(iPts.owner,point,iPts.points);
        }
        // notify
        PointNeighbors!(T) pn;
        pn.point=point;
        pn.neighbors=neighs;
        pn.dirDist=dirDist;
        pn.time=Clock.now;
        localContext.nCenter.notify("PointHasNewNeighbors",Variant(pn));
        return true;
    }

    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point.
    bool checkIfNeighbor(DynPVector!(T,XType)newPos,Point newPoint,T pSize) {
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
    /// point, add adds them to neighAtt and dirDist. If quickReturn is true returns immediately if
    /// it seems too far
    DirDistances!(T) addDirsOf(DynPVector!(T,XType)newPos,Point newPoint,T pSize,ref LocalGrowableArray!(PointAndDir) neighAtt,
        ref LocalGrowableArray!(DirDistances!(T)) dirDist)
    {
        bool veryClose=false;
        bool discarded=false;
        bool added=false;
        
        auto diff=newPos.dup();
        diff.opBypax(pos.dynVars.x,-1,1);
        DirDistances!(T) dDist;
        
        auto internalDiff=sqrt(diff.norm2()); // diff norm in internal coordinates
        dDist.xDist=internalDiff;
        if (dDist.veryFar(localContext)) return dDist;
        auto deriv1=pos.dynVars.dVarStruct.emptyDx();
        deriv1[]=0;
        pos.addToTSpace!(T)(diff,deriv1);
        diff.giveBack();
        // make an optimized version when the overlap is the identity??
        auto deriv1Dual=pos.dynVars.dVarStruct.emptyDualDx();
        pos.toDualTSpace!(T)(deriv1,deriv1Dual);
        // original norm in dual space
        auto origNorm=sqrt(deriv1Dual.norm2());
        dDist.dualDist=origNorm;
        if (dDist.veryFar(localContext)) return dDist;
        if (origNorm<localContext.minNormDual()*explorationSize){
            veryClose=true;
        }
        // norm in cartesian units (T space approx)
        auto normCartesian=sqrt(pos.dotInTSpace(deriv1,deriv1Dual));
        dDist.cartesianDist=normCartesian;
        if (dDist.veryFar(localContext)) return dDist;
        if (normCartesian<localContext.minRealNormSelf2){
            veryClose=true;
        }
        if (!dDist.neighbor(localContext,explorationSize,pSize)) return dDist;
        if (hasFrameOfRef){
            uint dirMax=invalidDir;
            size_t iMax;
            if (origNorm>localContext.minNormDual*explorationSize && origNorm<localContext.maxNormDual*explorationSize){ // this direction is the "main" direction and should be the first
                // it is a neighbor in some direction
                DynPVector!(T,DualDxType) deriv2Dual;
                /// checking which direction of the dual space we are really exploring
                synchronized(this){
                    if (minDir.isDummy){
                        deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
                    }
                }
                if (! deriv2Dual.isDummy()){
                    volatile T maxVal=0;
                    iMax=0;
                    foreach (i,v;deriv2Dual.pLoop){
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
                    bool noDir=false;
                    auto cosVal=maxVal/origNorm;
                    if (cosVal>localContext.sameDirCosAngle()){
                        auto dOpt2=pow2(explorationSize-maxVal)+pow2(origNorm)-pow2(maxVal);
                        if (dOpt2<localContext.dirDualSize2){
                            auto dOpt=((dOpt2>0)?cast(T)sqrt(dOpt2):cast(T)0);
                            DirDistances!(T) dDist1=dDist;
                            neighAtt(PointAndDir(newPoint,dirMax));
                            added=true;
                            dDist1.dualDirDist=dOpt;
                            dDist1.dualCos=cosVal;
                            dirDist(dDist1);
                            auto actualDirFlags=exploredDirs.atomicCAS(iMax,DirFlags.Explored,DirFlags.Free);
                        }
                    }
                    deriv2Dual.giveBack();
                }
            }
        
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
                lenDirs=repeat(ones!(T)(1),ndim)[0]; // should maybe give a function to construct this directly...
            }
            unaryOpStr!("*aPtr0=sqrt(*aPtr0);")(lenDirs);
            if (dirMax!=invalidDir){
                // complete data of dDist
                auto idealLen=explorationSize*lenDirs[iMax];
                auto dDist1=dirDist[dirDist.length-1];
                dDist1.cartesianCos=deriv2[iMax]/normCartesian;
                auto dOpt2=pow2(normCartesian)-pow2(deriv2[iMax])
                    +localContext.inDirCartesianScale2()*pow2(abs(deriv2[iMax])-idealLen);
                dDist1.cartesianDirDist=dOpt2;
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
                            exploredDirs.atomicCAS(dirAtt,DirFlags.Explored,DirFlags.Free);
                            if (idealLen<=localContext.zeroLen()){ // don't store directions for invalid dirs
                                DirDistances!(T) dDist1=dDist;
                                dDist1.cartesianCos=deriv2[idim]/normCartesian;
                                dDist1.cartesianDirDist=((dOpt2>0)?sqrt(dOpt2):0);
                                added=true;
                                neighAtt(PointAndDir(newPoint,dirAtt));
                                dirDist(dDist1);
                            }
                        }
                    }
                }
            }
        }
        if (veryClose || origNorm<=localContext.minNormDual*explorationSize){
            added=true;
            neighAtt(PointAndDir(newPoint,0));
            dirDist(dDist);
        }
        if ((!added)){
            if (origNorm<=localContext.maxNormDual*explorationSize+localContext.zeroLen()){
                added=true;
                neighAtt(PointAndDir(newPoint,ndirs)); // not in a specific direction
                dirDist(dDist);
            } else {
                assert(origNorm<=localContext.maxNormDual*pSize+localContext.zeroLen());
                added=true;
                neighAtt(PointAndDir(newPoint,invalidDir)); // a neighbor of the other point
                dirDist(dDist);
            }
        }
        return dDist;
    }
    /// adds the given point as neighbor
    void addNeighbor(Point p){
        synchronized(this){
            if (findFirstPred(neighbors.data,delegate bool(PointAndDir pDir){ return pDir.point==p; })!=neighbors.length){
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
        if (neighAtt.length==0){
            neighAtt(PointAndDir(p,invalidDir));
            dirDist(res);
        }
        addNeighbors(neighAtt.data,dirDist.data,hasGrad);
    }
    /// evaluates with the given context, returns if an evaluate was really done
    bool evalWithContext(LocalCalculationContext c){
        bool calcE=(gFlags&GFlags.EnergyEvaluated)==0;
        bool calcF=((gFlags&GFlags.GradientEvaluated)==0 && (localContext.cheapGrad()||(gFlags&GFlags.GradientInProgress)!=0));
        if (calcE||calcF){
            Real e=pos.dynVars.potentialEnergy;
            mixin(withPSys(`pSys[]=pos;`,"c."));
            c.updateEF(calcE,calcF);
            if (calcF) pos.checkMddx();
            mixin(withPSys(`pos[]=pSys;`,"c."));
            if (isNaN(pos.dynVars.potentialEnergy)) {
                pos.dynVars.potentialEnergy=e; // replace with if(!calcE) ?
            } else {
                e=pos.dynVars.potentialEnergy;
            }
            if (isNaN(e)){
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy after calculation in point ")(point);
                }),__FILE__,__LINE__);
            }
            auto target=localContext.ownerOfPoint(point);
            if (!calcF){
                localContext.addEnergyEvalLocal(target,point,e);
            }
            bool calcMore=false;
            if (!calcF){
                calcMore=localContext.speculativeGradient(target,point,e);
            } else {
                addEnergyEvalLocal(e);
            }
            if (calcMore){
                calcF=true;
                pos.checkMddx();
                c.updateEF(false,true);
                mixin(withPSys(`pos.dynVars.mddx[]=pSys.dynVars.mddx;`,"c."));
                // update also the energy? it should not change...
            }
            if (calcF){
                if ((gFlags&GFlags.LocalCopy)!=0){
                    localContext.addGradEvalLocal(target,point,pSysWriter(pos));
                } else {
                    didGradEval();
                }
            }
            return true;
        } else {
            return false;
        }
    }
    /// notifies a GFlags change from the given gFlags
    void notifyGFlagChange(uint oldGFlags){
        if (gFlags!=oldGFlags){
            GFlagsChange gFlagsChange;
            gFlagsChange.oldGFlags=oldGFlags;
            gFlagsChange.newGFlags=gFlags;
            gFlagsChange.point=point;
            if ((gFlags&GFlags.LocalCopy)==0){
                localContext.nCenter.notify("localPointChangedGFlags",Variant(gFlagsChange));
            } else {
                localContext.nCenter.notify("copiedPointChangedGFlags",Variant(gFlagsChange));
            }
        }
    }
    /// the energy of a neighbor was calculated
    void localNeighEnergy(PointAndDir[] neigh,Real e){
        typeof(_gFlags) oldGFlags;
        synchronized(this){
            if (isNaN(energy) || (gFlags&GFlags.EnergyEvaluated)==0) {
                assert((gFlags&GFlags.EnergyInfo)==GFlags.InProgress,collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy and not in progress in point ")(this); }));
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("skipping energy info (")(e)(") of neighbor ")(neigh)(" of point ")(point)
                        (" as this point doesn't have a calculated energy yet");
                });
                return;
            }
            oldGFlags=gFlags;
            if (e<energy && (gFlags&GFlags.NeighValsDecrease)==0){
                _gFlags|=GFlags.NeighValsDecrease;
            }else if (e>energy){
                _gFlags|=GFlags.NeighValsIncrease;
            } else {
                _gFlags|=GFlags.NeighValsSame;
            }
        }
        // check for special directions
        uint lastDir=ndirs-1;
        for(size_t iNeigh=0;iNeigh<neigh.length;++iNeigh){
            if (neigh[iNeigh].dir==1){ // along the force
                synchronized(this){
                    if (e<=energy) _gFlags|=GFlags.AlongForcesDecrease;
                    if (e>=energy) _gFlags|=GFlags.AlongForcesIncrease;
                }
            }
            if (neigh[iNeigh].dir==lastDir){ // along the gradient
                synchronized(this){
                    if (e<=energy) _gFlags|=GFlags.AlongGradientDecrease;
                    if (e>=energy) _gFlags|=GFlags.AlongGradientIncrease;
                }
            }
        }
        notifyGFlagChange(oldGFlags);
    }
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    void addEnergyEvalOther(Point p,Real e){
        PointAndDir[128] buf;
        LocalGrowableArray!(PointAndDir) localCopy=lGrowableArray(buf,0);
        // find neighbor
        synchronized(this){
            size_t iStart=0;
            while (1){
                while(iStart<neighbors.length){
                    if (p==neighbors[iStart].point) break;
                    ++iStart;
                }
                if (iStart==neighbors.length) break;
                size_t iEnd=iStart+1;
                while(iEnd<neighbors.length){
                    if (p!=neighbors[iEnd].point) break;
                    ++iEnd;
                }
                localCopy.appendArr(neighbors.data[iStart..iEnd]);
                iStart=iEnd;
            }
        }
        if (localCopy.length!=0){
            localNeighEnergy(localCopy.data,e);
        }
        localCopy.deallocData();
    }
    
    /// adds the energy to this point
    /// energy must be valid
    /// gFlags is updated at the end getting the energy of all neighbors
    void addEnergyEvalLocal(Real e){
        if (isNaN(e)){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("energy should not be NAN in point ")(point);
            }),__FILE__,__LINE__);
        }
        uint oldGFlags;
        synchronized(this){
            if ((! isNaN(pos.dynVars.potentialEnergy))&&(pos.dynVars.potentialEnergy!=e)){
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("double eval of energy (")(e)(" and ")(pos.dynVars.potentialEnergy)(") in point ")(point);
                });
            }
            pos.dynVars.potentialEnergy=e;
            oldGFlags=gFlags;
            _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
        }
        notifyGFlagChange(oldGFlags);
        updateNeighE();
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
        localNeigh.deallocData();
        auto energies=lMem.allocArr!(Real)(nPts.points.length);
        uint newGFlags=0;
        foreach (i;nPts.pLoop){
            auto silo=i.owner;
            localContext.neighborHasEnergy(silo,point,i.points,energy);
            auto energ=localContext.energyForPointsLocal(silo,i.points,energies[i.localLb..i.localUb]);
            auto nPts=i.points.length;
            for(size_t iPts=0;iPts<nPts;++iPts){
                addEnergyEvalOther(i.points()[iPts],energies[i.localLb+iPts]);
            }
        }
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags;
            _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyEvaluated;
        }
        notifyGFlagChange(oldGFlags);
    }
    
    /// the gradient of a neighbor was calculated (and possibly also the energy for the first time)
    void localNeighGrad(PointAndDir[] neigh,LazyMPLoader!(T)mainPoint,Real e){
        typeof(_gFlags) oldGFlags;
        synchronized(this){
            if (isNaN(energy) || (gFlags&GFlags.EnergyEvaluated)==0) {
                assert((gFlags&GFlags.EnergyInfo)==GFlags.InProgress,collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy and not in progress in point ")(this); }));
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("skipping energy info (")(e)(") of neighbor ")(neigh)(" of point ")(point)
                        (" as this point doesn't have a calculated energy yet");
                });
                return;
            }
            oldGFlags=gFlags;
            if (e<energy && (gFlags&GFlags.NeighValsDecrease)==0){
                _gFlags|=GFlags.NeighValsDecrease;
            }else if (e>energy){
                _gFlags|=GFlags.NeighValsIncrease;
            } else {
                _gFlags|=GFlags.NeighValsSame;
            }
        }
        // check for special directions
        uint lastDir=ndirs-1;
        for(size_t iNeigh=0;iNeigh<neigh.length;++iNeigh){
            if (neigh[iNeigh].dir==1){ // along the force
                synchronized(this){
                    if (e<=energy) _gFlags|=GFlags.AlongForcesDecrease;
                    if (e>=energy) _gFlags|=GFlags.AlongForcesIncrease;
                }
                if (hasFrameOfRef){
                    // to do: this is the place to put gradient based topology info, will be useful with free energy
                    //auto nPoint=mainPoint.mainPoint(localContext);
                    //nPoint.pos.dynVars.mddx // scalar prod with minDir, expected to be >0 (if it is a "normal" point), comparing with minDirScale gives an idea of the curvature
                }
            }
            if (neigh[iNeigh].dir==lastDir){ // along the gradient
                synchronized(this){
                    if (e<=energy) _gFlags|=GFlags.AlongGradientDecrease;
                    if (e>=energy) _gFlags|=GFlags.AlongGradientIncrease;
                }
                if (hasFrameOfRef){
                    // to do: this is the place to put gradient based topology info, will be useful with free energy
                    //auto nPoint=mainPoint.mainPoint(localContext);
                    //nPoint.pos.dynVars.mddx // scalar prod with minDir, expected to be >0 (if it is a "normal" point), comparing with minDirScale gives an idea of the curvature
                }
            }
        }
        notifyGFlagChange(oldGFlags);
    }
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    void addGradEvalOther(LazyMPLoader!(T)mp,Real e){
        PointAndDir[128] buf;
        LocalGrowableArray!(PointAndDir) localCopy=lGrowableArray(buf,0);
        Point p=mp.point;
        // find neighbor
        synchronized(this){
            size_t iStart=0;
            while (1){
                while(iStart<neighbors.length){
                    if (p==neighbors[iStart].point) break;
                    ++iStart;
                }
                if (iStart==neighbors.length) break;
                size_t iEnd=iStart+1;
                while(iEnd<neighbors.length){
                    if (p!=neighbors[iEnd].point) break;
                    ++iEnd;
                }
                localCopy.appendArr(neighbors.data[iStart..iEnd]);
                iStart=iEnd;
            }
        }
        if (localCopy.length!=0){
            localNeighGrad(localCopy.data,mp,e);
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
        bool newE=false;
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags;
            if ((gFlags&GFlags.EnergyEvaluated)==0){
                _gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
                newE=true;
            }
            _gFlags=(gFlags & ~GFlags.GradientInfo)|GFlags.GradientKnown;
        }
        buildMinDir();
        PointAndDir[128] buf;
        auto localNeigh=lGrowableArray(buf,0);
        synchronized(this){
            localNeigh.appendArr(neighbors.data);
            neighbors.clearData();
        }
        auto lData=localNeigh.data;
        auto nPts=new PointToOwner(&localContext.ownerOfPoint,delegate Point(size_t i){ return lData[i].point;},
            lData.length);
        Real[] energies;
        uint newGFlags=0;
        auto tNow=ev_time();
        auto thisP=LazyMPLoader!(T)(point,true,tNow);
        foreach (i;nPts.pLoop){
            auto silo=i.owner;
            localContext.neighborHasGradient(silo,thisP,i.points, energy);
            auto nPts=i.points.length;
            for(size_t iPts=0;iPts<nPts;++iPts){
                auto pointAtt=i.points()[iPts];
                scope(exit){ pAtt.release(); }
                auto pAtt=localContext.createLocalPoint(pointAtt,tNow);
                checkIfNeighbor(pAtt.pos.dynVars.x, pointAtt,pAtt.explorationSize);
                if (!isNaN(pAtt.energy)){
                    if (pAtt.hasFrameOfRef){
                        auto otherP=LazyMPLoader!(T)(pointAtt,true,tNow);
                        addGradEvalOther(otherP,pAtt.energy);
                        otherP.release();
                    } else {
                        addEnergyEvalOther(pointAtt,pAtt.energy);
                    }
                }
            }
        }
        notifyGFlagChange(oldGFlags);
        auto ownr=owner;
        localContext.notifyLocalObservers(delegate void(ExplorationObserverI!(T)obs){
            obs.addGradEvalLocal(ownr,point,pSysWriter(pos));
        });
    }
    /// builds the direction toward the minimum
    void buildMinDir(){
        assert((gFlags&GFlags.GradientInfo)==GFlags.GradientKnown,"needs gradient");
        auto newMinDir=pos.dynVars.dVarStruct.emptyDualDx();
        pos.toDualTSpace(pos.dynVars.mddx,newMinDir);
        auto n2=newMinDir.norm2();
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
        pos[]=p.pos;
        minDir[]=p.minDir;
        _exploredDirs=p.exploredDirs;
        _neighbors=p.neighbors;
        _minDirScale=p.minDirScale;
        _explorationSize=p.explorationSize;
        _gFlags=p.gFlags;
    }
}

MainPoint!(Real) dummyP;