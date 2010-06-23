/// this module represents 
module dchem.pnet.PointNetwork;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import Atomic=blip.sync.Atomic;
import blip.t.math.Math:max;
import tango.util.container.more.Heap;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import dchem.input.RootInput;
import blip.container.Deque;
import dchem.pnet.DirArray;
import dchem.pnet.PNetModels;

enum :uint{ invalidDir=uint.max; }

/// unique key for a point, and a direction in the same structure
struct PointAndDir{
    /// 16 bit direction, 36 bit sequence number, 12 bit collision reduction, the bit grouping might change, don't count on it
    ulong data; /// key of a main point (64 bit as trade off between collision probability and being small). 0 is an invalid key
    static PointAndDir opCall(Point k,uint dir){
        PointAndDir res;
        if (dir!=uint.max){
            assert(dir<0xF_FFFF,"direction too large");
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)dir)<<44);
        } else {
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
        return res;
    }
    /// direction
    uint dir(){
        uint res=cast(uint)(data>>44);
        return ((res!=0xF_FFFF)?res:uint.max);
    }
    //sets direction
    void dir(uint d){
        if (d!=uint.max){
            assert(d<0xF_FFFF,"direction too large");
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)d)<<44);
        } else {
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
    }
    /// base key
    Point point(){
        return Point(0x0000_0FFF_FFFF_FFFFUL & data);
    }
    /// sets base key, dropping extra bits
    void point(Point k){
        data=(0xFFFF_F000_0000_0000 & data)|(k.data & 0x0000_0FFF_FFFF_FFFFUL);
    }
    mixin(serializeSome("dchem.PointAndDir","data"));// split in Point and dir???
    mixin printOut!();
    struct ExpandedPointAndDir{
        PointAndDir pDir;
        uint dir(){
            return pDir.dir();
        }
        void dir(uint d){
            pDir.dir=d;
        }
        Point point(){
            return pDir.point();
        }
        void point(Point p){
            pDir.point=p;
        }
        mixin(serializeSome("dchem.PointDir","point|dir"));// split in Point and dir???
        mixin printOut!();
    }
    ExpandedPointAndDir expanded(){
        ExpandedPointAndDir e;
        e.pDir=*this;
        return e;
    }
}

/// a structure to keep point and its energy together
struct PointAndEnergy{
    Point point;
    Real energy;
    mixin(serializeSome("PointAndEnergy","point|energy"));
    mixin printOut!();
}

/// structure that collects togheter all points of the same owner (useful to bcast chuncks)
class PointToOwner{
    /// points sorted by owner & point
    Point[] points;
    /// index of the point in the original array
    size_t[] idx;
    /// list of owners
    SKey[] owner;
    /// start of the points of each owner (+ end, owner[i] has points[starts[i]..starts[i+1]])
    size_t starts;

    /// internal help structure (to reorder stuff)
    struct PIdxOwn{
        Point point;
        size_t idx;
        SKey owner;
    }
    /// returns a PointToOwner that contains the points extracted from pts ordered by owner and point
    /// if compressPts is true (default) then duplicate points are removed
    /// if keepIdx is true (the default) the indexes of the points in the origina array are kept
    this(T)(HashMap!(Points,SKey) ownerMap,T[]pts,bool compressPts=true,bool keepIdx=true){
        if (pts.length==0) return;

        PIdxOwn[128] buf;
        PIdxOwn[] points=buf;
        size_t i=0,j=1;
        points[0].point=pts[0].point;
        points[0].idx=0;
        points[0].owner=ownerMap[pts[0].point];
        while(j<pts.length){
            if ((!compressPts) || toSort[i].point!=pts[j].point){
                if (i==points.length){
                    auto newP=new PIdxOwn[](min(growLength(i+1,Point.sizeof),pts.length));
                    newP[0..i]=points;
                    if (i!=buf.length){
                        delete points;
                    }
                    points=newP;
                    idx=newIdx;
                }
                points[++i].point=pts[j].point;
                points[i].idx=j;
                points[i].owner=ownerMap[pts[j].point];
            }
            ++j;
        }
        ++i;
        points=points[0..i];
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
        if (i>buf.length){
            delete points.ptr;
        }
    }
    /// number of owners
    size_t nOwners(){
        return owner.length;
    }
    /// gets the i-th owner and its points
    bool ptsI(size_t i,out SKey owner,out Point[]points){
        if (i<owner.length){
            owner=this.owner[i];
            points=this.points[starts[i]..starts[i+1]];
            return true;
        }
        return false;
    }
    /// gets the i-th owner and its points and indexes
    bool ptsI(size_t i,out SKey owner,out Point[]points,out size_t idx){
        if (i<owner.length){
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
            return pOwn.points[starts[iOwn]..starts[iOwn+1]];
        }
        size_t[] idx(){
            return pOwn.idx[starts[iOwn]..starts[iOwn+1]];
        }
        size_t[] localLb(){
            return pOwn.starts[iOwn];
        }
        size_t[] localUb(){
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
        int opApply(int delegate(Iter*) loopBody){
            auto pool=cachedPoolNext(function Iter*(PoolI!(Iter*)p){
                auto res=new Iter;
                res.pool=p;
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
                for (size_t iOwn=0;iOwn<owner.length;++iOwn){
                    nIter=pool.getObj();
                    nIter.action=doIter;
                    nIter.iOwn=iOwn;
                    Task("PointToOwnerPLoopIter",&nIter.doAction).appendOnFinish(&nIter.giveBack).autorelease.submitYield();
                }
            }).autorelease.onFinish(&pool.stopCaching).executeNow();
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
struct PointNeighbors{
    Point point;
    PointAndDir[] neighbors;
    mixin(serializeSome("dchem.PointNeighbors","point|neighbors"));
    mixin printOut!();
}

/// distances of a point from a direction
struct DirDistances(T){
    T dualDist;
    T cartesianDist;
    T dualCos;
    T cartesianCos;
    T dualDirDist;
    T cartesianDirDist;
    /// all values are set
    bool full(){
        return !(isNAN(dualDist)||isNAN(cartesianDist)||
            isNAN(dualCos)||isNAN(cartesianCos)|| 
            isNAN(dualDirDist)||isNAN(cartesianDirDist));
    }
    bool veryFar(MainPoint p){
        if (!(isNAN(dualDist)||isNAN(cartesianDist)){
            return false;
        }
        if (dualDist>p.dualDistMax && cartesianDist>p.cartesianDistMax){
            return true;
        }
        return false;
    }
    bool neighbor(MainPoint p){
        if (!(isNAN(dualDist)||isNAN(cartesianDist)){
            return false;
        }
        if (dualDist>p.dualDistMax && cartesianDist>p.cartesianDistMax){
            return false;
        }
        return true;
    }
}

/// a main evaluation point
class MainPoint(T){
    /// encodes a discrete probability level
    enum Prob{
        Unlikely=0, /// from the current information it is unlikely that the point has that feature
        Possible=1, /// one indicator is pointing toward the point having that feature
        Probable=2, /// it is likely that the point has that feature
        Confirmed=3, /// it has been confirmed that the point has that feature
    }
    /// flags of the point
    enum GFlags{
        None=0,/// no flags
        EnergyInfo =(3<<0),
        EnergyInProgress=(1<<0),        /// energy calculation at this point is in progress
        EnergyKnown =(3<<0),              /// the energy is known
        EnergyEvaluated =(2<<0),        /// the energy was evaluated (b-cast done)
        GradientInfo =(3<<2),
        GradientInProgress=(1<<2),      /// gradient calculation at this point is in progress
        GradientKnown =(3<<2),            /// gradient is known
        GradientEvaluated =(2<<2),      /// gradient was evaluated, frame of reference is established, b-cast done
        InProgress=EnergyInProgress|GradientInProgress; /// a calculation is in progress
        DoNotExplore=(1<<4),            /// no further exploration should be performed starting from this point
        FullyExplored=(1<<5),           /// all directions around this point have been explored
        FullyEvaluated=(1<<6),          /// all directions around this point have been explored and evaluated
        OldApprox=(1<<7),               /// energy/gradient are based on old data, newer better data should be available
        LocalCopy=(1<<8),               /// this MainPoint is a local copy done for efficency reasons (localContext is not the owner)
        AlongForces=(3<<9),             /// codify what happens going along the minimuma direction
        AlongForcesUnknow=0,            /// not yet claculated
        AlongForcesDecrease=(1<<9),     /// along the forces the energy decreases (normal point)
        AlongForcesIncrease=(2<<9),     /// along the forces the energy increases (critical point close)
        AlongGradient=(3<<11),          /// codifies what happens along the gradient
        AlongGradientUnknow=0,          /// not yet calculated
        AlongGradientDecrease=(1<<11),  /// along the gradient the energy decreases (critical point close)
        AlongGradientIncrease=(2<<11),  /// along the gradient the energy increases (normal point)
        NeighVals=(7<<13),              /// codifies the energy values of the neighbors
        NeighValsUnknown=0,             /// no neighbor values known yet
        NeighValsDecrease=(1<<13),      /// some neighbors are smaller
        NeighValsSame=(2<<13),          /// some neighbors have the same energy
        NeighValsIncrease=(4<<13),      /// some neighbors are larger
        AttractorBorder=(3<<16),        /// codifies the current attractor status
        AttractorBorderForces=(1<<16),  /// we are at the border going along the forces
        AttractorBorderGradient=(2<<16),/// we are at the border going along the gradients
    }
    LocalSilosI!(T) localContext;
    /// identification of this point
    Point point;
    /// attractor of this point (attractors are identified by their minima)
    Point attractor;
    /// position of the point (energy, derivatives,... might be invalid, only real positions have to be valid)
    ParticleSys!(T) pos;
    /// direction of the minimum in the dual space with norm 1 wrt. euclidean norm
    /// (frame of reference for minimization), valid only if this point is a starting point for further exploration
    DynPVect!(T,2) minDir;
    /// bit array of the directions that have been explored (is of length null all direction have been explored)
    /// this stores only real exploration directions, directions that a given method blocks should be stored in
    /// a separated DirArray
    DirArray exploredDirs;
    /// neighbors, and the direction with respect to this point in which they are
    LocalGrowableArray!(PointAndDir) neighbors;
    /// neighbor distances, stored only if requested. 6 numbers for each neighbor:
    /// dualDist, cartDist, cosDual, rDistDual, cartCos, cartRDist
    LocalGrowableArray!(T) neighDistances;
    /// scale of mindir to recover the dual gradient (useful??)
    T minDirScale;
    /// exploration size for this point (used to establish neighbors,...) 
    T explorationSize;
    /// bit-or of GFlags of the current point
    uint gFlags;
    // *** topology, functions of gFlags
    /// probability of having a minimum close by for the given gFlags
    static Prob minimumForGFlags(uint gFlags){
        if ((gFlags&GFlags.NeighValsDecrease)!=0) return Prob.Unlikely;
        /// the next one accepts all elements of a flat potential as minima, which while logially correct overcrowds the special points
        /// add requirement that AlongForcesDecrease is 0 ???
        if ((gFlags&GFlags.AlongForcesIncrease)!=0) return (((gFlags&GFlags.FullyEvaluated)!=0)?Prob.Confirmed:Prob.Probable);
        if ((gFlags&GFlags.NeighValsIncrease)!=0) return Prob.Possible;
        return Prob.Unlikely;
    }
    /// probability of having a critical point close by for the given gFlags
    static Prob criticalPointForGFlags(uint gFlags){
        if ((gFlags&GFlags.AlongForcesIncrease)!=0) return Prob.Probable;
        if ((gFlags&GFlags.AlongGradientDecrease)!=0) return Prob.Probable;
        return Prob.Unlikely;
    }
    /// probability of having a saddle point
    Prob saddlePointForGFlags(uint gFlags){
        if (criticalPointForGFlags(gFlags)==Prob.Probable){
            if ((gFlags&GFlags.AttractorBorder)!=0){
                return Prob.Probable;
            }
            if ((gFlags&GFlags.NeighValsIncrease)!=0 && (gFlags&GFlags.NeighValsDecrease)!=0){
                return Prob.Likely;
            }
            return Prob.Possible;
        }
        return Prob.Unlikely;
    }
    /// probability that close to this point there is a special (critical) point
    Prob specialPointForGFlags(uint gFlags){
        int p1=cast(int)minimumForGFlags(gFlags), p2=cast(int)criticalPointForGFlags(gFlags),
            p3=cast(int)saddlePointForGFlags(gFlags);
        return cast(Prob)(max(p1,max(p2,p3)));
    }
    /// different point types, considered as exclusive types
    enum PointType{
        Minimum,
        TransitionPoint,
        TransitionSurface,
        CriticalPoint,
        NormalPoint,
    }
    /// returns the main that this point can be considered as
    PointType pointTypeForGFlags(uint gFlags){
        if (minimumForGFlags(gFlags)>Prob.Possible){
            return PointType.Minimum;
        } else if (saddlePointForGFlags(gFlags)>Prob.Possible){
            return PointType.TransitionPoint;
        } else if (criticalPointForGFlags(gFlags)>Prob.Possible){
            return PointType.CriticalPoint;
        } else if ((gFlags&GFlags.AttractorBorder)!=0){
            return PointType.TransitionSurface;
        } else {
            return PointType.NormalPoint;
        }
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
    
    /// structure to encode a main point efficienly (useful to transfer it o another context)
    struct DriedPoint{
        Point point;
        PSysWriter!(T) pos;
        DynPVectorWriter!(T,DualDxType) minDir;
        DirArray exploredDirs;
        LocalGrowableArray!(PointAndDir) neighbors;
        T minDirScale;
        T explorationSize;
        uint gFlags;
        
        static DriedPoint opCall(MainPoint p){
            DriedPoint res;
            res.point=p.point;
            res.pos=pSysWriter(p.pos);
            res.minDir=dynPVectorWriter(p.minDir);
            res.exploredDirs=p.exploredDirs;
            res.neighbors=p.neighbors;
            res.minDirScale=minDirScale;
            res.explorationSize=explorationSize;
            return res;
        }
        
        mixin(serializeSome("dchem.DriedPoint","point|pos|minDir|exploredDirs|neighbors|minDirScale|explorationSize|gFlags"));
        mixin printOut!();
        MainPoint hydrate(LocalSilosI!(T)localContext){
            auto pSys=localContext.refPos.dup();
            pSys[]=pos;
            DynPVect!(T,2) mDir;
            if (minDir.isNonNull()){
                mDir=minDir.toPVector(pSys.dualDxGroup,true);
            }
            MainPoint res=new MainPoint(localContext,point,pSys,gFlags|GFlags.LocalCopy,explorationSize,
                mDir,minDirScale,exploredDirs,neighbors);
            return res;
        }
    }
    /// this point as dried point
    DriedPoint driedPoint(){
        return DriedPoint(this);
    }
    
    /// constructor
    this(MinEExplorer!(T) localContext,Point point,ParticleSys!(T) pos,uint gFlags=GFlags.None,
        T explorationSize=-1,DynPVect!(T,2) minDir=null,T minDirScale=T.init,
        FlagsArray exploredDirs=null,PointAndDir[] neighbors=null)
    {
        this.localContext=localContext;
        this.point=point;
        this.pos=pos;
        this.gFlags=gFlags;
        this.explorationSize=explorationSize;
        if (this.repulsionSize<0) this.repulsionSize=localContext.repulsionSize;
        if (this.explorationSize<=0) this.explorationSize=localContext.explorationSize;
        this.minDir=minDir;
        this.minDirScale=minDirScale;
        this.exploredDirs=exploredDirs;
        if (this.exploredDirs is null){
            this.exploredDirs=new FlagsArray(nDirs);
        }
        this.neighbors=lGrowableArray(neighbors,neighbors.length,GASharing.Global);
    }
    
    /// explores the next direction, immediately marks it as in exploration
    /// returns direction 0 if the gradient has to be evaluated, an invalid direction if the current point is in evaluation,
    /// and an invalid point only if all direction are explored/ the point should not be explored
    /// if lastIsLast is true then the last direction (gradient) is explored as last
    PointAndDir exploreNext(DirFlags methodDirs=null,bool lastIsLast=true){
        auto exploreFlags=DirFlags.Explored;
        if ((gFlags&GFlags.DoNotExplore)!=0){
            return PointAndDir(0,uint.max);
        }
        if ((gFlags&GFlags.GradientEvaluated)==0){
            synchronized(this){
                if ((gFlags&(GFlags.InProgress|GFlags.GradientEvaluated))==0){
                    gFlags|=GFlags.GradientInProgress;
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
                    size_t nextDir=dirFlags.length;
                    if(dirFlags.findFreeAndSet(1,exploreFlags)){
                        nextDir=0;
                    } else {
                        auto nextDir=exploredDirs.findFreeAndSet(start,exploreFlags,);
                    }
                    if (nextDir!=dirFlags.length) return PointAndDir(this.point,nextDir,lastIsLast,methodDirs);
                    gFlags|=GFlags.FullyExplored;
                }
            }
        }
        return PointAndDir(0,0);
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
    DynPVector!(T,DualDxType) dualDir(uint dir){
        if (dir==0) return pos;
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        T val=(neg?-1:1);
        auto v0=minDir.emptyCopy;
        v0[]=0;
        v0[rDir]=val;
        auto newDir=rotateEiV(V,M)(0,minDir,v0);
        return newDir;
    }
    /// energy of the current point
    Real energy(){
        synchronized(this){
            return pos.dynVars.potentialEnergy;
        }
    }
    /// calculates the position exploring from here in the given direction
    /// returns null
    ParticleSys!(T) posInDir(uint dir){
        auto newDir=dualDir(dir);
        auto cartesianNorm=pos.projectInTSpace(newDir);
        auto projNorm=newDir.norm2();
        if (cartesianNorm<=localContext.minRealNormSelf0){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)("discarded because cartesian norm is too small. Projected norm was ")(projNorm)(".");
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
        pos.addFromDualTSpace(T)(newDir,resPSys.dynVars.x);
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
    /// and returns a valid Point for it.
    bool acceptNewDirection(Point newPoint,uint dir){
        // Point
        if (dir==0){
            throw new Exception("should not call acceptNewDirection when exploring the point itself (direction 0)",
                __FILE__,__LINE__);
        }
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        auto localPoint=localContext.localPo
        /// verify pos:
        auto diff=newPos.dup();
        diff.axpby(pos.dynVars.x,-1,1);
        if (localContext.minMoveInternal>0 && sqrt(diff.norm2())<localContext.minMoveInternal){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in the internal coordinates of exploration direction ")(dir)(" of point ")
                    (point)(" is too small in internal coordinates ")(sqrt(diff.norm2()))(", discarding evaluation");
            });
            return false;
        }
        auto deriv1=pos.dynVars.emptyDx();
        deriv1[]=0;
        pos.addToTSpace(diff,deriv1);
        diff.giveBack();
        // make an optimized version when the overlap is the identity??
        auto deriv1Dual=pos.toDualTSpace(deriv1,deriv1Dual);
        // original norm in dual space
        auto origNorm=sqrt(deriv1Dual.norm2());
        if (origNorm<=localContext.minNormDualSelf()*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is too small:")(origNorm)(", discarding evaluation");
            });
            return false;
        }
        if (origNorm>maxNormDual*explorationSize || origNorm<localContext.minNormDual*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is outside the exploration zone with length ")(origNorm)("*")(explorationSize)
                    (", continuing evaluation");
            });
        }
        // norm in cartesian units
        normCartesian=sqrt(pos.dotInTSpace(deriv1,deriv1Dual));
        if (normCartesian<localContext.minRealNormSelf2){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in real space of exploration in direction ")(dir)(" of point ")(point)
                    (" is too small:")(normCartesian)(", discarding evaluation");
            });
        }
        // go in the minDir frame of reference
        auto deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
        auto deriv2=rotateVEi(minDir,0,deriv1);
        
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);

        checkDir(deriv2,deriv2Dual,normCartesian,normDual,neighAtt);
        
        if (find(neighAtt.data,delegate bool(PointAndDir p){ return p.dir==dir; })==neighAtt.length){
            // not in the expected direction
            if (neighAtt.length==0){ // not a neighbor...
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                    (" is not a neighbor, continuing evaluation..."); // stop instead?
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
                
                auto dirMax=toDir(iMax,deriv2Dual[iMax]<0);
                bool nonVisited=false;
                for (size_t idir=0;idir<neighAtt.length;++idir){
                    if (neighAtt[idir].dir==0){
                        if (dirMax==dir){
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is  mostly in the expected direction in the dual space, but ")
                                    (iMax)(", but in the expected region, continuing evaluation");
                            });
                            nonVisited=true;
                        } else {
                            localContext.logMsg(delegate void(CharSink s){
                                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                                    (" is  very close to the starting point, and not in the expected direction");
                        }
                    } else {
                        auto actualDirFlags=exploredDirs.atomicCAS(iMax,GFlags.Explored,GFlags.Free);
                        if (actualDirFlags == GFlags.Free){
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
        deriv2Dual.giveBack();
        deriv2.giveBack();
        if (neighbors.length!=0){
            synchronized(this){
                assert(find(neighbors.data,delegate bool(PointAndDir p){ return p.point==newPoint; })==neighbors.length,
                    "point "~p.toString()~" was already added to neighbors of "~point.toString());
                neighbors.appendArr(neighAtt.data); // duplicates are stored contiguously
            }
            
            point.addNeighbor
        neighAtt.deallocData();
        return newPoint;
    }
    
    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point.
    bool checkNeighbors(DynPVector!(T,XType)newPos,Point newPoint,T pSize) {
        // ignore if the point is already in the neighbors
        synchronized(this){
            if (find(neighbors.data,delegate bool(PointAndDir p){ return p.point==newPoint; })!=neighbors.length){
                /// already added
                return true;
            }
        }
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);
        bool veryClose=false;
        bool discarded=false;
        
        auto diff=newPos.dup();
        diff.axpby(pos.dynVars.x,-1,1);
        auto internalDiff=sqrt(diff.norm2()); // diff norm in internal coordinates
        if (internalDiff<localContext.minMoveInternal){
            veryClose=true;
        } else if (localContext.maxMoveInternal>0 && internalDiff>localContext.maxMoveInternal) {
            // too far discarding
            discarded=true;
            version(Quick){
                return false;
            }
        }
        auto deriv1=pos.dynVars.emptyDx();
        deriv1[]=0;
        pos.addToTSpace(diff,deriv1);
        diff.giveBack();
        // make an optimized version when the overlap is the identity??
        auto deriv1Dual=pos.toDualTSpace(deriv1,deriv1Dual);
        // original norm in dual space
        auto origNorm=sqrt(deriv1Dual.norm2());
        if (origNorm<localContext.minNormDual()*explorationSize){
            veryClose=true;
        }
        // norm in cartesian units (T space approx)
        normCartesian=sqrt(pos.dotInTSpace(deriv1,deriv1Dual));
        if (normCartesian<localContext.minRealNorm){
            veryClose=true;
        }
        if (origNorm>minNormDual*explorationSize && origNorm<maxNormDual*explorationSize){ // this direction is the "main" direction and should be the first
            // it is a neighbor in some direction
            DynPVector!(T,DualDxType) deriv2Dual;
            /// checking which direction of the dual space we are really exploring
            synchronized(this){
                if (minDir is null){
                    deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
                }
            }
            if (! deriv2Dual.isDummy()){
                volatile T maxVal=0;
                size_t iMax;
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
                auto dirMax=toDir(iMax,deriv2Dual[iMax]<0);
                bool noDir=false;
                auto cosVal=maxVal/origNorm;
                if (cosVal>localContext.sameDirCosAngle()){
                    auto dOpt2=pow(localContext.optimalNormDual*explorationSize-maxVal,2)+pow(origNorm,2)-pow(maxVal,2);
                    auto dOpt=((dOpt2>0)?cast(T)sqrt(dOpt2):cast(T)0);
                    if (dOpt<dirDualSize){
                        neighAtt.append(PointAndDir(newPoint,dirMax));
                        auto actualDirFlags=exploredDirs.atomicCAS(iMax,GFlags.Explored,GFlags.Free);
                    }
                }
                deriv2Dual.giveBack();
            }
        }
        
        // check for ill defined directions in cartesian space
        auto deriv2=rotateVEi(minDir,0,deriv1);
        // length of the directions in cartesian units
        auto ov=pos.maybeOverlap();
        scope rotOv=rotateVEi(minDir,0,ov.data.T.dup(true));
        auto lenDirs=diag(rotOv);
        unaryOpStr!("*aPtr0=sqrt(*aPtr0);")(lenDirs);
        // loop on all directions checking cartesian thresholds
        auto dim=ndim;
        for(idim=0;idim<dim;++idim){
            auto idealLen=explorationSize*lenDirs[i];
            if (abs(normCartesian-idealLen)<minRealNorm &&
                abs(deriv2[idim])>=localContext.sameDirCosAngleCartesian()*explorationSize*lenDirs[i])
            {
                auto dOpt2=pow(normCartesian,2)-pow(deriv2[idim],2)+pow(deriv2[idim]-idealLen,2)
                if (dOpt2=<dirCartesianSize2){
                    auto dirAtt=toDir(idim,deriv2[idim]<0);
                    dirFlags.atomicCAS(dirAtt,DirFlags.Explored,DirFlags.Free);
                    if (deriv2[idim]!=0){ // don't store directions for invalid dirs
                        neighAtt.append(PointAndDir(newPoint,dirAtt));
                    }
                }
            }
        }
        
        if (veryClose || origNorm<=minNormDual*explorationSize){
            neighAtt.approx(PointAndDir(newPoint,0));
        } else if (origNorm<=maxNormDual*explorationSize && neighAtt.length==0){
                neighAtt.append(PointAndDir(newPoint,nDirs)); // not in a specific direction
            }
        }
        if (origNorm<localContext.maxNormDual*pSize){
            neighAtt.append(PointAndDir(newPoint,invalidDir)); // a neighbor of the other point
        } else {
            return false;
        }
        synchronized(this){
            /// find if we already added this point...
            auto pos=find(neighbors.data,delegate bool(PointAndDir pD){ return pD.point==newPoint; });
            if (pos!=neighbors.data.length){
                bool same=true;
                auto nLength=neighbors.data.lengthM
                if (neighAtt.length+pos>=neighbors.length) {
                    same=false;
                } else {
                    for (size_t ipos=0;ipos<neighAtt.length;++ipos){
                        if (neighbors.data[pos+ipos]!=neighAtt.data[ipos]){
                            same=false;
                            break;
                        }
                    }
                }
                if (same){
                    localContext.logMsg(delegate void(CharSink s){
                        dumper(s)("double add of neighbor ")(newPoint)(" to ")(point)(", ignored");
                    })
                } else {
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        size_t lastPos=pos;
                        while(lastPos<neighbors.length && neighbors[lastPos]==neighbors[pos]){
                            ++lastPos;
                        }
                        dumper(s)("double add of neighbor ")(newPoint)(" to ")(point)(" generated different neighbors:")
                            (neighAtt.data)(" vs ")(neighbors.data[pos..lastPos]);
                    }),__FILE__,__LINE__);
                }
            }
            neighbors.appendArr(neighAtt.data); // duplicates are stored contiguously
        }
        // communicate to neighbors
        auto neigToNotify=new PointToOwner(neighAtt.data);
        foreach (iPts;neigToNotify.pLoop){
            localContext.silos[iPts.owner].areNeigh(point,iPts.points);
        }
        // notify
        PointNeighbors pn;
        pn.point=point;
        pn.neighbors=neighAtt.takeData();
        localContext.nCenter.notify("PointHasNewNeighbors",Variant(pn));
        assert(!discarded,"did discard valid point");
        return true;
    }
    /// adds the given point as neighbor
    void addNeighbor(Point p){
        synchronized(this){
            if (find(neighbors.data,delegate bool(PointAndDir p){ return p.point==newPoint; })!=neighbors.length){
                return;
            }
        }
        auto op=localContext.ownerOfPoint(p);
        DynPVector!(T,XType) newPos=localContext.silos[op].pointPosLocal(p).dynPVector(pos)
        if(!checkNeighbors(newPos,Point newPoint,this.explorationSize)){// to do what do we really need to pass? what should be done about cartesian thresholds?
            synchronized(this){
                neighbors.append(PointAndDir(p,invalidDir));
            }
        }
    }
    /// evaluates with the given context, returns if an evaluate was really done
    bool evalWithContext(CalculationContext c){
        bool calcE=(gFlags&GFlags.EnergyEvaluated)==0;
        bool calcF=((gFlags&GFlags.GradientEvaluated)==0 && (cheapGrad()||gFlags&GFlags.GradientInProgress));
        if (calcE||calcF){
            Real e=pSys.dynVars.potentialEnergy;
            mixin(withPSys(`pSys[]=pos;`,"c."));
            c.updateEF(calcE,calcF);
            mixin(withPSys(`pos[]=pSys;`,"c."));
            if (isNAN(pos.dynVars.potentialEnergy)) {
                pos.dynVars.potentialEnergy=e; // replace with if(!calcE) ?
            } else {
                e=pos.dynVars.potentialEnergy;
            }
            auto target=localContext.activeExplorers[localContext.owner];
            if (isNAN(e)){
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("invalid energy after calculation in point ")(point);
                }),__FILE__,__LINE__);
            }
            bool calcMore=false;
            if (!calcF){
                calcMore=target.speculativeGradient(point,e);
            }
            if (calcMore){
                calcF=true;
                c.updateEF(false,true);
                if (isNAN(pos.dynVars.potentialEnergy)) {
                    pos.dynVars.potentialEnergy=e; // replace with if(!calcE) ?
                } else {
                    e=pos.dynVars.potentialEnergy;
                }
            }
            if (calcF){
                if ((gFlags&GFlags.LocalCopy)!=0){
                    assert(target is activeExplorers[localContext.key],"non local copy and not local");
                    target.addGradEvalLocal(point,pSysWriter(pos));
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
        typeof(gFlags) oldGFlags;
        synchronized(this){
            if (isNAN(energy)) {
                assert((gFlags&GFlags.InProgress)!=0,"invalid energy and not in progress");
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("skipping energy info (")(e)(") of neighbor ")(neighbor)(" of point ")(point)
                        (" as this point doesn't have a calculated energy yet");
                });
                return;
            }
            auto oldGFlags=gFlags;
            if (e<energy && (gFlags&GFlags.NeighValsDecrease)==0){
                gFlags|=GFlags.NeighValsDecrease;
            }else if (e>energy){
                gFlags|=GFlags.NeighValsIncrease;
            } else {
                gFlags|=GFlags.NeighValsSame;
            }
        }
        // check for special directions
        uint lastDir=ndirs-1;
        for(iNeigh=0;iNeigh<neigh.length;++iNeigh){
            if (neigh[iNeigh].dir==1){ // along the force
                synchronized(this){
                    if (e<=energy) gFlags|=GFlags.AlongForcesDecrease;
                    if (e>=energy) gFlags|=GFlags.AlongForcesIncrease;
                }
            }
            if (neigh[iNeigh]==lastDir){ // along the gradient
                synchronized(this){
                    if (e<=energy) gFlags|=GFlags.AlongGradientDecrease;
                    if (e>=energy) gFlags|=GFlags.AlongGradientIncrease;
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
                localCopy.appendArr(neighbors[iStart..iEnd]);
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
        if (isNAN(e)){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("energy should not be NAN in point ")(point);
            }),__FILE__,__LINE__);
        }
        uint oldGFlags;
        synchronized(this){
            if ((! isNAN(pos.dynVars.potentialEnergy))&&(pos.dynVars.potentialEnergy!=e)){
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("double eval of energy (")(e)(" and ")(pos.dynVars.potentialEnergy)(") in point ")(point);
                });
            }
            pos.dynVars.potentialEnergy=e;
            oldGFlags=gFlags;
            gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
        }
        notifyGFlagChange(oldGFlags);
        updateNeighE();
    }
    /// gets the energies of the neighbors
    void updateNeighE(){
        PointAndDir[128] buf;
        Real[128] buf2;
        auto lMem=LocalMem(buf2);
        auto localNeigh=lGrowableArray(buf,0);
        synchronized(this){
            localNeigh.appendArr(neighbors.data);
        }
        auto nPts=PointToOwner(localNeigh.data);
        localNeigh.deallocData();
        auto energies=lMem!(Real)(nPts.points.length);
        uint newGFlags=0;
        foreach (i,nPts.pLoop){
            auto energ=silos[owner].energiesOfPoints(i.points,energies[i.localLb,i.localUb]);
            auto nPts=i.points.length;
            for(size_t iPts=0;iPts<nPts;++iPts){
                addEnergyEvalOther(i.points()[iPts],energies[i.localLb+iPts]);
            }
        }
        uint oldGFlags;
        synchronized(this){
            oldGFlags=gFlags
            gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyEvaluated;
        }
        notifyGFlagChange(oldGFlags);
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
                gFlags=(gFlags & ~GFlags.EnergyInfo)|GFlags.EnergyKnown;
                newE=true;
            }
            gFlags=(gFlags & ~GFlags.GradientInfo)|GFlags.GradientKnown;
        }
        buildMinDir();
        PointAndDir[128] buf;
        Real[128] buf2;
        auto lMem=LocalMem(buf2);
        auto localNeigh=lGrowableArray(buf,0);
        synchronized(this){
            localNeigh.appendArr(neighbors.data);
            neighbors.clearData();
        }
        auto nPts=PointToOwner(localNeigh.data);
        Real[] energies;
        if (newE) energies=lMem!(Real)(nPts.points.length);
        uint newGFlags=0;
        foreach (i,nPts.pLoop){
            Real[] energ;
            energ=silos[owner].energiesOfPoints(i.points,energies[i.localLb,i.localUb]);
            auto nPts=i.points.length;
            for(size_t iPts=0;iPts<nPts;++iPts){
                auto pointAtt=i.points()[iPts];
                auto pPos=silos[owner].pointPos();
                auto pPos2=pPos.toDynPVector(localContext.refPos.dynVars.dVarStruct.xGroup);
                checkNeighbors(pPos2, pointAtt);
                addEnergyEvalOther(pointAtt,energ[iPts]);
            }
        }
        notifyGFlagChange(oldGFlags);
    }
    /// builds the direction toward the minimum
    void buildMinDir(){
        assert((gFlags&GFlags.GradientEvaluated),"needs gradient");
        auto newMinDir=pos.dynVars.emptyDualDx();
        pos.toDualTSpace(pos.dx,newMinDir);
        auto n2=newMinDir.norm2();
        if (n2==0){
            newMinDir[0]=1;
            minDirScale=0
        } else {
            minDirScale=sqrt(n2);
            newMinDir*=cast(T)(1/minDirScale);
        }
        synchronized(this){
            minDir=newMinDir;
        }
    }
    
    void retain(){
        assert(refCount!=0,"refCount was 0 in retain");
        atomicAdd(refCount,1);
    }
    void release(){
        assert(refCount!=0,"refCount was 0 in retain");
        atomicAdd(refCount,-1);
        if (refCount==0){
            if ((gFlags&GFlags.LocalCopy)!=0){
                
            }
        }
    }
}

