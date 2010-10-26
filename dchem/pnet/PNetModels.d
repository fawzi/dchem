/// interfaces of a PointNetwork
module dchem.pnet.PNetModels;
import blip.serialization.Serialization;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.sys.DynVars;
import dchem.sys.Constraints;
import dchem.pnet.DirArray;
import dchem.input.RootInput;
import dchem.Common;
import dchem.input.WriteOut;
import blip.parallel.smp.Wait;
import blip.util.RefCount;
import blip.container.GrowableArray;
import blip.time.Time;
import blip.time.Clock;
import blip.util.NotificationCenter;
import dchem.calculator.CalculatorModels;
import blip.math.IEEE;
import blip.math.random.Random;
import blip.container.Cache;
import blip.container.Pool;

/// unique key for a point, dimensions might change
/// key might have some special values,
///   0: means no point,
///   1: pseudo key for points whose exploration has been excluded (for example due to constraints)
/// the key itself has the first 10 bits (size might change) that are used for collision reduction
/// and the rest that encodes a mostly sequential number (this is indirectly useful as the exploration is not random),
/// and to allow efficient storage of dense arrays (making unbundling of properties of points easy)
/// all keys with 0 at the beginning and only the first 10 bit set are thus never generated,
/// and might be used as special keys (see the special values above)
struct Point{
    ulong data;
    ulong key(){
        return data;
    }
    /// sets base key (should not be too large)
    void key(ulong k){
        assert(k<=0x0000_0FFF_FFFF_FFFFUL,"key too large");
        data=k&0x0000_0FFF_FFFF_FFFFUL; // drop the extra &?
    }
    /// set base key without loosing too many bits (and guranteeing that if it was non 0 it stays non 0)
    static Point opCall(ulong k){
        Point res;
        res.data=(k & 0x0000_0FFF_FFFF_FFFFUL);
        auto r=((k>>21)&0x0000_07FF_FF80_0000UL);
        if (r!=0){
            res.data=(res.data^r)|0x0000_0800_0000_0000UL;
        }
        return res;
    }
    bool isValid(){
        return (data&0x0000_0FFF_FFFF_F000UL)!=0;
    }
    mixin(serializeSome("dchem.Point","data"));
    mixin printOut!();
    hash_t toHash(){
        static if(hash_t.sizeof==4){
            return cast(hash_t)(0xFFFF_FFFF&(data^(data>>32)));
        } else {
            return cast(hash_t)data;
        }
    }
    int opCmp(Point p2){
        return ((data>p2.data)?1:((data==p2.data)?0:-1));
    }
}

enum :uint{ invalidDir=uint.max }

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
    int opCmp(PointAndDir p2){
        return ((data>p2.data)?1:((data==p2.data)?0:-1));
    }
}

/// a structure to keep point and its energy together
struct PointAndEnergy{
    Point point;
    Real energy;
    mixin(serializeSome("PointAndEnergy","point|energy"));
    mixin printOut!();
    int opCmp(PointAndEnergy p2){
        return ((energy<p2.energy)?-1:((energy>p2.energy)?1:point.opCmp(p2.point)));
    }
}

/// structure to encode a main point efficienly (useful to transfer it o another context)
struct DriedPoint(T){
    Point point;
    PSysWriter!(T) pos;
    DynPVectorWriter!(T,DualDxType) minDir;
    FlagsArray exploredDirs;
    LocalGrowableArray!(PointAndDir) neighbors;
    T minDirScale;
    T explorationSize;
    uint gFlags;
    
mixin(serializeSome("dchem.DriedPoint!("~T.stringof~")","point|pos|minDir|exploredDirs|neighbors|minDirScale|explorationSize|gFlags"));
    mixin printOut!();
}

/// distances of a point from a direction
struct DirDistances(T){
    T xDist; /// internal (x space) distance
    T dualDist; /// distance in the dual space
    T cartesianDist; /// cartesianDistance (in the TS)
    T dualCos; /// cosinus with the direction in the dual space
    T cartesianCos; /// cosinus with the direction in the cartesian metric
    T dualDirDist; /// distance in the dual space from the perfect direction point
    T cartesianDirDist; /// distance in the cartesian metric from the perfect direction point
    /// all values are set
    bool full(){
        return !(isNaN(dualDist)||isNaN(cartesianDist)||
            isNaN(dualCos)||isNaN(cartesianCos)|| 
            isNaN(dualDirDist)||isNaN(cartesianDirDist));
    }
    bool veryFar(LocalSilosI!(T) silos){
        if (((!isNaN(xDist)) && silos.maxMoveInternal()>0 && xDist>silos.maxMoveInternal())||
            ((!isNaN(dualDist)) && silos.maxMoveDual()>0 && dualDist>silos.maxMoveDual())||
            ((!isNaN(cartesianDist)) && silos.maxMoveCartesian()>0 && cartesianDist>silos.maxMoveCartesian()))
        {
            return true;
        }
        return false;
    }
    bool neighbor(LocalSilosI!(T) silos,Real eDist1,Real eDist2){
        if (veryFar(silos)) return false;
        if (((!isNaN(dualDist)) && (dualDist<=eDist1*silos.maxNormDual || dualDist<=eDist2*silos.maxNormDual))||
            ((!isNaN(cartesianDist)) && cartesianDist<=silos.maxNormCartesian)){
            return true;
        }
        return true;
    }
    mixin(serializeSome("dchem.DirDistances!("~T.stringof~")",`
        xDist: internal (x space) distance
        dualDist: distance in the dual space
        cartesianDist: cartesianDistance (in the TS)
        dualCos: cosinus with the direction in the dual space
        cartesianCos: cosinus with the direction in the cartesian metric
        dualDirDist: distance in the dual space from the perfect direction point
        cartesianDirDist: distance in the cartesian metric from the perfect direction point
    `));
    mixin printOut!();
}

/// encodes a discrete probability level
enum Prob{
    Unlikely=0, /// from the current information it is unlikely that the point has that feature
    Possible=1, /// one indicator is pointing toward the point having that feature
    Likely=2, /// it is likely that the point has that feature
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
    InProgress=EnergyInProgress|GradientInProgress, /// a calculation is in progress
    DoNotExplore=(1<<4),            /// no further exploration should be performed starting from this point
    FullyExplored=(1<<5),           /// all directions around this point have been explored
    FullyEvaluated=(1<<6),          /// all directions around this point have been explored and evaluated
    OldApprox=(1<<7),               /// energy/gradient are based on old data, newer better data should be available (also set when dropped)
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
    PointBcasted=(1<<17),           /// the point was broadcasted around (dropping is a global op)
}

// *** topology, functions of gFlags
/// probability of having a minimum close by for the given gFlags
Prob minimumForGFlags(uint gFlags){
    if ((gFlags&GFlags.NeighValsDecrease)!=0) return Prob.Unlikely;
    /// the next one accepts all elements of a flat potential as minima, which while logially correct overcrowds the special points
    /// add requirement that AlongForcesDecrease is 0 ???
    if ((gFlags&GFlags.AlongForcesIncrease)!=0) return (((gFlags&GFlags.FullyEvaluated)!=0)?Prob.Confirmed:Prob.Likely);
    if ((gFlags&GFlags.NeighValsIncrease)!=0) return Prob.Possible;
    return Prob.Unlikely;
}
/// probability of having a critical point close by for the given gFlags
Prob criticalPointForGFlags(uint gFlags){
    if ((gFlags&GFlags.AlongForcesIncrease)!=0) return Prob.Likely;
    if ((gFlags&GFlags.AlongGradientDecrease)!=0) return Prob.Likely;
    return Prob.Unlikely;
}
/// probability of having a saddle point
Prob saddlePointForGFlags(uint gFlags){
    if (criticalPointForGFlags(gFlags)==Prob.Likely){
        if ((gFlags&GFlags.AttractorBorder)!=0){
            return Prob.Likely;
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
/// returns the main PointType that this point can be considered as
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

/// a main evaluation point
interface MainPointI(T){
    /// local context (silos) in which the point is stored
    LocalSilosI!(T) localContext();
    /// identification of this point
    Point point();
    /// attractor of this point (attractors are identified by their minima)
    Point attractor();
    /// position of the point (energy, derivatives,... might be invalid, only real positions have to be valid)
    ParticleSys!(T) pos();
    /// direction of the minimum in the dual space with norm 1 wrt. euclidean norm
    /// (frame of reference for minimization), valid only if this point is a starting point for further exploration
    DynPVector!(T,DualDxType) minDir();
    /// bit array of the directions that have been explored (is of length null all direction have been explored)
    /// this stores only real exploration directions, directions that a given method blocks should be stored in
    /// a separated DirArray
    FlagsArray exploredDirs();
    /// neighbors, and the direction with respect to this point in which they are
    LocalGrowableArray!(PointAndDir) neighbors();
    /// neighbor distances, stored only if requested. 6 numbers for each neighbor:
    /// dualDist, cartDist, cosDual, rDistDual, cartCos, cartRDist
    LocalGrowableArray!(DirDistances!(T)) neighDistances();
    /// scale of mindir to recover the dual gradient (useful??)
    T minDirScale();
    /// exploration size for this point (used to establish neighbors,...) 
    T explorationSize();
    /// bit-or of GFlags of the current point
    uint gFlags();
    /// atomic cas on the flags of this point
    uint gFlagsAtomicCAS(uint newVal,uint oldVal);
    /// atomic op on the flags of this point
    uint gFlagsAtomicOp(uint delegate(uint) op);
    /// if this point is a local copy, and not the "main" point
    bool isLocalCopy();
    /// if the point is explorable
    bool isExplorable();
    /// if the point is still valid (ie. not superseded or dropped)
    bool isValid();
    /// it the point was broadcasted to all silos
    bool isPublic();
    /// if the point has a local frame of reference (i.e. minDir)
    bool hasFrameOfRef();
    
    // *** topology, utility methods
    /// probability of having a minimum close by (utility method)
    Prob minimum();
    /// probability of having a critical point close by (utility method)
    Prob criticalPoint();
    /// probability of having a saddle point (utility method)
    Prob saddlePoint();
    /// probability that close to this point there is a special (critical) point (utility method)
    Prob specialPoint();
    /// returns the main that this point can be considered as (utility method)
    PointType pointType();
    
    /// returns the number of dimensions
    uint ndim();
    /// returns the number of directions (counting also 0, the core dir)
    uint ndirs();
    /// returns the dir value for the given dimension and sign
    uint toDir(uint idim,bool neg);
    /// transforms a non core (i.e. non 0) dir to dimension and sign
    void fromDir(uint dir,out uint dim,out bool neg);
    
    /// this point as dried point
    DriedPoint!(T) driedPoint();
    
    /// explores the next direction, immediately marks it as in exploration
    /// returns direction 0 if the gradient has to be evaluated, an invalid direction if the current point is in evaluation,
    /// and an invalid point only if all direction are explored/ the point should not be explored
    /// if lastIsLast is true then the last direction (gradient) is explored as last
    PointAndDir exploreNext(FlagsArray methodDirs=null,bool lastIsLast=true);
    /// a point in the given direction in the dual space
    DynPVector!(T,DualDxType) createDualDir(uint dir);
    /// energy of the current point
    Real energy();
    /// calculates the position exploring from here in the given direction
    /// returns null if the position is *really* too close
    ParticleSys!(T) createPosInDir(uint dir);
    /// returns if the given direction is acceptable. If it is accepted sets the local directions covered
    /// & neighbors, and broadcasts the point.
    /// this is most likely called on a copy of the point, and where newPoint is local...
    bool acceptNewDirection(Point newPoint,uint dir);
    
    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point.
    bool checkIfNeighbor(DynPVector!(T,XType)newPos,Point newPoint,T pSize) ;
    /// adds the given point as neighbor
    void addNeighbor(Point p);
    /// internal method, adds the pre-processed points&dirs to the list of neighbors and notifies neighbors
    /// returns true if the points were really added
    bool addNeighbors(PointAndDir[] neighs,DirDistances!(T)[] dirDists,bool hadGrad);
    
    /// evaluates with the given context, returns if an evaluation was really done
    bool evalWithContext(CalculationContext c);
    /// notifies a GFlags change from the given gFlags
    void notifyGFlagChange(uint oldGFlags);
    /// the energy of a neighbor was calculated
    void localNeighEnergy(PointAndDir[] neigh,Real e);
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    void addEnergyEvalOther(Point p,Real e);
    /// the gradient of a neighbor was calculated (and possibly also the energy for the first time)
    void localNeighGrad(PointAndDir[] neigh,LazyMPLoader!(T)mainPoint,Real e);
    /// adds an energy evaluation for the given point (that should be a neighbor of this point)
    void addGradEvalOther(LazyMPLoader!(T)p,Real e);
    
    /// adds the energy to this point
    /// energy must be valid
    /// gFlags is updated at the end getting the energy of all neighbors
    void addEnergyEvalLocal(Real e);
    /// gets the energies of the neighbors and broadcasts its own energy to them
    void updateNeighE();
    /// to be called after the gradient was added to the pos of this point
    /// energy must be valid
    /// creates dualDir, minDirScale, updates gFlags, assign directions to all neighbors
    /// might modify the attractor, and notifies neighbors
    void didGradEval();
    /// marks this point as dropped, returns true if the point was not yet dropped
    bool drop();
    
    void retain();
    void release();
    /// assign from a dired point
    void opSliceAssign(DriedPoint!(T));
}

/// lazy loader of points
class LazyMPLoader(T):Serializable{
    Point point;
    MainPointI!(T) _mainPoint;
    WaitCondition waitPoint;
    Time maxTime;
    bool weakUpdate;
    Exception exception;
    mixin RefCountMixin!();
    PoolI!(LazyMPLoader) pool;
    void release0(){
        if (_mainPoint!is null) _mainPoint.release();
        _mainPoint=null;
        if (pool!is null) pool.giveBack(this);
    }
    bool hasPoint(){
        return _mainPoint!is null;
    }
    // if needed loads the mainPoint, otherwise returns the cached value
    MainPointI!(T)mainPoint(LocalSilosI!(T) silos){
        if(_mainPoint is null){
            if (exception){
                throw exception;
            }
            bool shouldLoad=false;
            synchronized(this){
                if (waitPoint is null){
                    waitPoint=new WaitCondition(&hasPoint);
                    shouldLoad=true;
                }
            }
            if (shouldLoad){
                auto mp=silos.createLocalPoint(point,maxTime);
                assert(_mainPoint is null);
                try{
                    _mainPoint=mp;
                } catch(Exception e){
                    exception=e;
                }
                waitPoint.notifyAll();
                if (exception!is null) throw exception;
                return _mainPoint;
            } else {
                waitPoint.wait();
            }
            if (exception!is null) throw exception;
            assert(_mainPoint!is null);
        }
        if (weakUpdate){
            auto mp=silos.createLocalPoint(point,maxTime);
            synchronized(this){
                if (mp!is _mainPoint){
                    auto oldP=_mainPoint;
                    _mainPoint=mp;
                    oldP.release();
                } else {
                    mp.release();
                }
                return _mainPoint;
            }
        }
        return _mainPoint;
    }
    /// creates a LazyMPLoader that will load the given point from silos
    /// if weakUpdate is true tries to get the newest version that is locally available
    /// the time is used to decide if the cached version is new enough, or a new version has to be
    /// fetched from the owning silos
    this(Point point,bool weakUpdate,Time time,PoolI!(LazyMPLoader) pool=null){
        this.point=point;
        this.weakUpdate=weakUpdate;
        this.maxTime=time;
        this.pool=pool;
    }
    /// ditto
    this(Point point,bool weakUpdate=true){
        this(point,weakUpdate,Clock.now);
    }
    this(){}
    typeof(this) postUnserialize(Unserializer s){
        maxTime=Clock.now; // reset time as it cannot be assumed to be synchronized between different computers
        return this;
    }
    mixin(serializeSome("dchem.LazyMPLoader!("~T.stringof~")",`point|weakUpdate|maxTime`));
    mixin printOut!();
    
    static PoolI!(LazyMPLoader) gPool;
    static this(){
        gPool=cachedPool(function LazyMPLoader(PoolI!(LazyMPLoader)p){
            return new LazyMPLoader(Point(0),true,Time(0),p);
        });
    }
    static LazyMPLoader opCall(Point point,bool weakUpdate,Time time,PoolI!(LazyMPLoader) pool=null){
        auto res=gPool.getObj();
        res.maxTime=time;
        res.weakUpdate=weakUpdate;
        res.point=point;
        return res;
    }
    static LazyMPLoader opCall(Point point,bool weakUpdate=true){
        return LazyMPLoader.opCall(point,weakUpdate,Clock.now);
    }
}

alias ulong SKey; /// silos key
enum SKeyVal:SKey{
    Invalid=0,
    Any=1,
    All=2
}

/// calls that make an exploration progress, and correspond to notifications that are called.
///
/// All methods begin with an SKey that is the "target" of the message, and could also be used
/// for redundancy/load balancing in the future. This approach make it easy to hide a complex
/// parallel network behind a single connection.
interface ExplorationObserverI(T){
    /// stops a silos
    /// at the moment there is no support for dynamic adding/removal of silos, (worth adding???)
    /// so s should be only SKeyVal.All
    void shutdown(SKey s);
    /// adds energy for a point local to s and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy);
    /// adds gradient value to a point that should be owned by s. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys);
    /// communicates to s that the given point is being explored
    /// pSize is the point size, flags the flags of the point
    void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags);
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to s
    void neighborHasEnergy(SKey s,Point p,Point[] neighbors,Real energy);
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, Real energy);
    
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point,SKey owner);
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source);
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey,Point);
}

/// an object that can offer new points to explore
/// actually should not inherit from ExplorationObserverI, but this way we avoid multiple inheritance bugs
/// these method are *not* public/remote
interface ExplorerI(T):ExplorationObserverI!(T){
    /// returns a point to evaluate, calls available when a new point might be available, returns 1 if one should wait, 0 if the explorer has finished exploring
    Point pointToEvaluateLocal(void delegate(ExplorerI!(T)) availableAgain);
    /// called when an evaluation fails, returns 1, retry, 0 no choice, -1 do not retry
    int evaluationFailedLocal(Point);
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradientLocal(Point p,Real energy);
}

/// a task that can be transferred to a silos and performed there
interface RemoteSilosOpI(T):Serializable{ // could be Serializable and SilosWorkerI, avoiding it due to compiler bugs
    void workOn(LocalSilosI!(T)silos);
    void stop();
}

/// interface of a silos (storage) of the point network
///
/// just like ExplorationObserverI all methods have SKey as first argument (see there for the rationale)
interface PNetSilosI(T):ExplorationObserverI!(T){
    /// returns a point to evaluate, returns 0 if there are no further points, 1 if there are no points now
    /// but there should be in the future
    Point pointToEvaluate(SKey);
    /// called when an evaluation fails
    void evaluationFailed(SKey s,Point);
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradient(SKey s,Point p,Real energy);

    /// load (usage) of the silos s in some arbitrary units
    Real load(SKey s);
    
    /// energy for the local points (NAN if not yet known)
    Real[] energyForPointsLocal(SKey s,Point[],Real[]);
    /// energy for the points (NAN if not yet known)
    Real[] energyForPoints(SKey s,Point[],Real[]);
    /// returns a snapshot of the given point (asking first silos s)
    DriedPoint!(T)mainPoint(SKey s,Point);
    /// returns a snapshot of the given point that is local to the silos s
    DriedPoint!(T)mainPointLocal(SKey s,Point);
    /// owner of the given point (asking s first)
    SKey pointOwner(SKey s,Point);
    /// the next free silos (for storage)
    SKey nextFreeSilos(SKey s);

    /// tells the local points neighs that they have p0 as neighbor
    void addPointToLocalNeighs(SKey s,Point p0,Point[]neighs);
    /// tells the local point p0 that neighDirs (computed when p0 gradient was hadGrad) should be added to it
    bool addNeighDirsToLocalPoint(SKey s,Point p0,PointAndDir[]neighDirs,DirDistances!(T)[]dirDists,bool hadGrad);
    /// operation to be executed on the given silos
    void executeLocally(SKey s,RemoteSilosOpI!(T) op);
    
    /// dictionary with the values of the various properties
    Real[char[]] propertiesDict(SKey s);
    // expose creation & bcast of points and update from dried points? merging should be done carefully to avoid problems... so for now you should do them via executeLocal...
}

const char[] silosMethodsStr=`shutdown|addEnergyEvalLocal|addGradEvalLocal|publishPoint|neighborHasEnergy|neighborHasGradient|finishedExploringPoint|didLocalPublish|publishCollision|pointToEvaluate|evaluationFailed|speculativeGradient|load|energyForPointsLocal|energyForPoints|mainPoint|mainPointLocal|pointOwner|nextFreeSilos|addPointToLocalNeighs|addNeighDirsToLocalPoint|executeLocally|propertiesDict`;

/// local interface, to a silos (basically a silos client)
interface LocalSilosI(T): PNetSilosI!(T) {
    /// if this silos is owner of the key k (note that using this rather than having a single key per silos
    /// will allow more complex interpretations of SKey in the future...)
    bool hasKey(SKey k);
    /// random number generator
    RandomSync rand();
    /// writes out a log message (at once)
    void logMsg(void delegate(void delegate(char[]))writer);
    /// writes out a log message
    void logMsg(char[]msg);
    /// owner of the given point (just a utility method)
    SKey ownerOfPoint(Point);
    
    /// adds an extra observer that will be notified about the network changes
    void addObserver(ExplorationObserverI!(T) o);
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o);
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o);
    /// removes the given explorer
    void rmExplorer(ExplorerI!(T)o);
    /// notify observers, the operation should not raise (or the whole program stops)
    void notifyLocalObservers(void delegate(ExplorationObserverI!(T))n);
    
    /// reference position, this should be used just to create other ParticleSystem, never directly
    ParticleSys!(T) refPos();
    /// constraints for this system
    ConstraintI!(T) constraints();
    /// if the gradient is cheap to compute
    bool cheapGrad();
    /// local notification center
    NotificationCenter nCenter();
    
    // values to define the point net topology
    
    // quick cutoffs
    /// if bigger than zero quickly discards all points that are further than this value in internal units
    T maxMoveInternal();
    /// if bigger than zero quickly discards all points that are further than this value in the cartesian space
    T maxMoveCartesian();
    /// if bigger than zero quickly discards all points that are further than this value in the dual space
    T maxMoveDual();
    
    // neighs
    /// points that are less than maxNormCartesian apart are considered neighbors
    T maxNormCartesian();
    /// points that are less than maxNormDual apart (in units of explorationSize) are cosidered neighbors
    /// (this is used also to accept exploratios)
    T maxNormDual();
    
    /// minimum move in internal units to accept a direction
    T minMoveInternal();
    /// minimum norm of uint vector after projection of null direction to be considered as a valid direction
    T minProjectionResidual();
    /// maximum cosinus of the angle between two vectors that are still considered in the same direction in dual space 
    T sameDirCosAngle();
    /// radius in which to look for neighbors (in units of explorationSize)
    T sameDirCosAngleCartesian();
    /// minimum norm in the dual T space to accept an exploration (in units of explorationSize)
    T minNormDual();
    /// minimum norm in the dual T space to accept an exploration for a self generated direction 
    /// (in units of explorationSize)
    T minNormDualSelf();
    /// minimum norm in the real (cartesian) space for a self generated direction to be accepted before moving
    T minRealNormSelf0();
    /// minimum norm in the real (cartesian) space to which a self generated direction should be rescaled
    T minRealNormSelf1();
    /// minimum norm in the real (cartesian) space for the real movement (after constraints,...)
    /// for a self generated direction to be accepted
    T minRealNormSelf2();
    /// explorationSize
    T explorationSize();
    /// square of the maximum distance from the optimal direction point in explorationSize units
    T dirDualSize2();
    /// special scale used in the cartesian space to rescale the direction perfectly in the direction direction
    T inDirCartesianScale2();
    /// square of maximum distance from the optimal direction point in cartesian units
    /// If the distance is x along the direction and y along the direction dirCartesianSize2 is compared
    /// to inDirCartesianScale2*x*x+y*y
    T dirCartesianSize2();
    /// length that is considered 0: the invalid directions should be orthogonal to the dualDir, so that
    /// the lenght along them should be exactly 0, but nomerical error might make them slightly larger than 0
    T zeroLen();
    
    /// creates a new point located at newPos in this silos, the point is not yet broadcasted
    /// not all silos might support creation of local points, use nextFreeSilos to get a silos
    /// thas supports it, use executeLocally to create, setup & publish a point...
    MainPointI!(T) newPointAt(DynPVector!(T,XType) newPos);
    /// local point mainpoint (the real reference point)
    MainPointI!(T) mainPointL(Point);
    /// publishes the local point given.
    /// Returns the pubblished point which might be different fron the argument if a collision did happen.
    MainPointI!(T) bcastPoint(MainPointI!(T));
    /// a local point (possibly a copy), is retained, and needs to be released (thus the create in the name)
    /// the time t is used to decide if a cached version can be used
    MainPointI!(T)createLocalPoint(Point p,Time t);
    /// drops a cached point (the point is not in use anymore)
    void dropCachedPoint(MainPointI!(T)p);

    /// an url that can be used to contact the silos core
    char[] silosCoreUrl();
}

/// an object that works on a LocalSilos
interface SilosWorkerI(T){
    void workOn(LocalSilosI!(T) silos);
}

/// a loader of points (what is created by the input)
interface SilosWorkerGen:InputElement{
    SilosWorkerI!(Real) silosWorkerReal();
    SilosWorkerI!(LowP) silosWorkerLowP();
}

/// define silosWorkerT extractor helper
mixin(genTypeTMixin("SilosWorker","silosWorker","",""));

/// an exploration observer generator (what is created by the input)
interface ExplorationObserverGen:InputElement{
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos);
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos);
}

/// define observerT extractor helper
mixin(genTypeTMixin("ExplorationObserver","observer","LocalSilosI!(T)silos","silos"));

/// an exploration generator (what is created by the input)
interface ExplorerGen:ExplorationObserverGen{
    ExplorerI!(Real) explorerReal(LocalSilosI!(Real) silos);
    ExplorerI!(LowP) explorerLowP(LocalSilosI!(LowP) silos);
}

/// define explorerT extractor helper
mixin(genTypeTMixin("Explorer","explorer","LocalSilosI!(T)silos","silos"));

/// list of the properties exposed by LocalSilosI (useful for ctfe)
const char[] propertiesList=`discretizationStep|minProjectionResidual|sameDirCosAngle|minNormDual|minNormDualSelf
    minRealNormSelf0|minRealNormSelf1|maxNormDual|explorationSize|dirSize2|dirCartesianSize2
    maxMoveInternal|maxMoveCartesian|maxMoveDual|maxNormCartesian|minMoveInternal|sameDirCosAngleCartesian|
    minRealNormSelf2|dirDualSize2|inDirCartesianScale2|zeroLen`;