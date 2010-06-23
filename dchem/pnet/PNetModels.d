/// interfaces of a PointNetwork
module dchem.pnet.PNetModels;
import blip.serialization.Serialization;
import blip.

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
}

/// calls that make an exploration progress, and correspond to notifications that are called
interface ExplorationObserverI(T){
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(Point p,Real energy);
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(Point p,PSysWriter!(T) pSys);
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void addExploredPoint(EKey owner,Point point,PSysWriter!(T) pos,uint flags);
    
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(Point);
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(Point);
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void dropPoint(EKey,Point);
}

/// an object that can offer new points to explore
/// actually should not inherit from ExplorationObserverI, but this way we avoid multiple inheritance bugs
interface ExplorerI(T):ExplorationObserverI!(T){
    /// returns a point to evaluate
    Point pointToEvaluate();
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(Point,uint flags);
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradient(Point p,Real energy);
}

alias ulong SKey; /// silos key

/// interface of a silos (storage) of the point network
interface PNetSilosI(T):ExplorationObserverI!(T){
    /// stops all silos (at the moment there is no support for dynamic adding/removal of silos, worth adding???)
    void shutdown();
    /// key of this silos
    SKey key();
    /// where the journal is kept
    char[] journalPos();
    /// load (usage) of the silos in some arbitrary units
    Real load();
    
    
    /// returns the position of the given local point
    PSysWriter!(T)pointPosLocal(Point);
    /// returns the position of the given point
    PSysWriter!(T)pointPos(Point);
    /// energy for the local points (NAN if not yet known)
    Real[] energyForPointsLocal(Point[],Real[]);
    /// energy for the points (NAN if not yet known)
    Real[] energyForPoints(Point[],Real[]);
    /// returns a snapshot of the given point
    DriedMainPoint!(T)mainPoint(Point);
    /// returns a snapshot of the given point
    DriedMainPoint!(T)mainPointLocal(Point);
    /// tells that the given points are neighbors of p0
    void addNeighLocal(Point p0,Point[]neighs);
    /// informs that source has processed point p0
    void processedLocal(Point p0,SKey source);
}

interface LocalSilosI!(T):PNetSilosI{
    /// adds an extra observer that will be notified about the network changes
    void addObserver(ExplorationObserverI!(T) o);
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o);
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o);
    /// removes the given explorer
    void rmExplorer(ExplorerI!(T)o);
    
    /// reference position, this should be used just to create other ParticleSystem, never directly
    ParticleSystem!(T) refPos;
    /// constraints for this system
    MultiConstraints constraints();
    /// if the gradient is cheap to compute
    bool cheapGrad();
    /// local notification center
    NotificationCenter nCenter;
    /// owner of the given point
    SKey ownerOfPoint(Point);
    /// local point mainpoint (the real reference point)
    MainPoint!(T) mainPointL(Point);
    /// creates a new point somwhere in the network
    Point newPointAt(DynPVector!(T,XType) newPos);
    
    // values to define the point net topology
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
    /// max norm in the dual T space to accept an exploration (in units of explorationSize)
    T maxNormDual();
    /// explorationSize
    T explorationSize();
    /// square of the maximum distance from the optimal direction point in explorationSize units
    T dirSize2();
    /// square of maximum distance from the optimal direction point in cartesian units
    T dirCartesianSize2();
    /// a local point (possibly a copy), is retained, and needs to be released
    MainPoint!(T)localPoint(Point p);
}