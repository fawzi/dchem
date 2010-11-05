/// a silos client
module dchem.pnet.PNetSilosClient;
import dchem.pnet.PNetModels;
import blip.parallel.rpc.RpcBase;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import dchem.input.RootInput;
import dchem.calculator.Calculator;
import dchem.calculator.CalculatorModels;
import dchem.pnet.MainPoint;
import blip.math.random.Random;
import dchem.sys.Constraints;
import blip.util.NotificationCenter;
import dchem.input.WriteOut;
import blip.math.IEEE;
import dchem.sys.DynVars;
import blip.parallel.smp.Wait;
import blip.io.Console;
import blip.parallel.smp.WorkManager;
import blip.container.Deque;
import blip.container.Cache;
import blip.container.Pool;

/// parameters for PNetSilosClient
class PNetSilosClientGen:InputElement {
    char[] connectionUrl;
    long ownerCacheSize=10;
    mixin(serializeSome("dchem.PNetSilosClient",`
    connectionUrl: the url to use to get the connection to the silos
    ownerCacheSize: amount of points whose owner is cached (10)`));
    mixin printOut!();
    mixin myFieldMixin!();
    
    this(char[] connectionUrl,long ownerCacheSize=10){
        this.connectionUrl=connectionUrl;
        this.ownerCacheSize=ownerCacheSize;
    }
    bool verify(CharSink s){
        bool res=true;
        if (ownerCacheSize<0){
            dumper(s)("ownerCacheSize in field ")(myFieldName)(" should be non negative\n");
            res=false;
        }
        return res;
    }
}

struct PointAndOwner{
    Point point;
    SKey owner;
}

char[] extractFromProperties(char[]props){
    char[] res;
    auto propsP=ctfeSplit("| \n",props,true);
    foreach(property;propsP){
        if (property.length>0){
            res~=`
    T `~property~`(){
        return cast(T)properties["`~property~`"];
    }`;
        }
    }
    return res;
}

class PNetSilosClient(T){
    PNetSilosClientGen input;
    PNetSilos!(T) connection;
    Real[char[]] properties;
    RandomSync _rand;
    CharSink log;
    Deque!(PointAndOwner) ownerCache;
    ParticleSys!(T) _refPos;
    ConstraintI!(T) _constraints;
    HashMap!(Point,CachedPoint!(T)) localCache;
    NotificationCenter _nCenter;
    CalculationContext ctx;
    
    this(PNetSilosClientGen input, CalculationContext ctx, PNetSilosI!(T) connection=null,NotificationCenter nCenter=null){
        this.input=input;
        this.ctx=ctx;
        if (ctx is null){
            assert(0,"CalculationContext should not be null");
            //input.evaluator.contentT!(Method)().getCalculator(true,[]);
        }
        if (_refPos is null){
            CalculationContext cInstance=ctx;
            switch(cInstance.activePrecision()){
            case Precision.Real:
                _refPos=cInstance.pSysReal.dupT!(T)(PSDupLevel.All);
                auto c1=cInstance.constraintsReal();
                if (c1 !is null)
                    _constraints=constraintT!(T)(c1.constraintGen,_refPos);
                break;
            case Precision.LowP:
                _refPos=cInstance.pSysLowP.dupT!(T)(PSDupLevel.All);
                auto c2=cInstance.constraintsLowP();
                if (c2 !is null)
                    _constraints=constraintT!(T)(c2.constraintGen,_refPos);
                break;
            default:
                assert(0,"unknown activePrecision");
            }
        }
        if (_constraints is null) _constraints=new NoConstraint!(T)();
        
        this.connection=connection;
        if (connection is null){
            this.connection=ProtocolHandler.proxyForUrlT!(PNetSilosI!(T))(input.connectionUrl);
        }
        this._nCenter=nCenter;
        if (nCenter is null){
            this._nCenter=new NotificationCenter();
        }
        _rand=new RandomSync();
        this.log=sout;
        properties=this.propertiesDict(SKeyVal.Any);
    }
    // ExplorationObserverI
    
    /// stops a silos
    /// at the moment there is no support for dynamic adding/removal of silos, (worth adding???)
    /// so s should be only SKeyVal.All
    void shutdown(SKey s,int speed){
        connection.shutdown(s,speed);
    }
    /// adds energy for a point local to s and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy){
        connection.addEnergyEvalLocal(s,p,energy);
    }
    /// adds gradient value to a point that should be owned by s. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){
        connection.addGradEvalLocal(s,p,pSys);
    }
    /// communicates to s that the given point is being explored
    /// pSize is the point size, flags the flags of the point
    void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
        connection.publishPoint(s,owner,point,pos,pSize,flags);
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to s
    void neighborHasEnergy(SKey s,Point p,Point[] neighbors,Real energy){
        connection.neighborHasEnergy(s,p,neighbors,energy);
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, Real energy){
        connection.neighborHasGradient(s,p,neighbors,energy);
    }
    
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point p,SKey owner){
        connection.finishedExploringPoint(s,p,owner);
    }
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source){
        connection.didLocalPublish(s,p0,source);
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey k,Point p){
        connection.publishCollision(k,p);
    }

    // an object that can offer new points to explore
    // actually should not inherit from ExplorationObserverI, but this way we avoid multiple inheritance bugs
    // these method are not public/remote
    // ExplorerI(T)
    
    /// returns a point to evaluate
    Point pointToEvaluate(SKey s){
        return connection.pointToEvaluate(s);
    }
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(SKey s,Point p){
        connection.evaluationFailed(s,p);
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradient(SKey s,Point p,Real energy){
        return connection.speculativeGradient(s,p,energy);
    }

    // interface of a silos (storage) of the point network
    //
    // just like ExplorationObserverI all methods have SKey as first argument (see there for the rationale)
    // PNetSilosI(T)
    
    /// load (usage) of the silos s in some arbitrary units
    Real load(SKey s){
        return connection.load(s);
    }
    
    /// energy for the local points (NAN if not yet known)
    Real[] energyForPointsLocal(SKey s,Point[] pts,Real[] ens){
        return connection.energyForPointsLocal(s,pts,ens);
    }
    /// energy for the points (NAN if not yet known)
    Real[] energyForPoints(SKey s,Point[] pts,Real[] ens){
        return connection.energyForPoints(s,pts,ens);
    }
    /// returns a snapshot of the given point (asking first silos s)
    DriedPoint!(T)mainPoint(SKey s,Point p){
        return connection.mainPoint(s,p);
    }
    /// returns a snapshot of the given point that is local to the silos s
    DriedPoint!(T)mainPointLocal(SKey s,Point p){
        return connection.mainPointLocal(s,p);
    }
    /// owner of the given point (asking s first, which can be useful for non public points)
    SKey pointOwner(SKey s,Point p){
        return connection.pointOwner(s,p);
    }
    /// the next free silos (for storage)
    SKey nextFreeSilos(SKey s){
        return connection.nextFreeSilos(s);
    }

    /// tells the local points neighs that they have p0 as neighbor
    void addPointToLocalNeighs(SKey s,Point p0,Point[]neighs){
        connection.addPointToLocalNeighs(s,p0,neighs);
    }
    /// tells the local point p0 that neighDirs (computed when p0 gradient was hadGrad) should be added to it
    bool addNeighDirsToLocalPoint(SKey s,Point p0,PointAndDir[]neighDirs,DirDistances!(T)[]dirDists,bool hadGrad){
        return connection.addNeighDirsToLocalPoint(s,p0,neighDirs,dirDists,hadGrad);
    }
    /// operation to be executed on the given silos
    void executeLocally(SKey s,RemoteSilosOpI!(T) op){
        connection.executeLocally(s,op);
    }
    
    /// dictionary with the values of the various properties
    Real[char[]] propertiesDict(SKey s){
        return connection.propertiesDict(s);
    }

    // local interface, to a silos (basically a silos client)
    // LocalSilosI(T)

    /// if this silos is owner of the key k (note that using this rather than having a single key per silos
    /// will allow more complex interpretations of SKey in the future...)
    bool hasKey(SKey k){
        return false;
    }
    /// random number generator
    RandomSync rand(){
        return _rand;
    }
    /// writes out a log message (at once)
    void logMsg(void delegate(void delegate(char[]))writer){
        sinkTogether(log,writer);
    }
    /// writes out a log message
    void logMsg(char[]msg){
        log(msg);
    }
    /// owner of the given point (just a utility method)
    SKey ownerOfPoint(Point p){
        PointAndOwner res;
        synchronized(ownerCache){
            ownerCache.filterInPlace(delegate bool(PointAndOwner po){
                if (po.point==p){
                    res=po;
                    return false;
                }
                return true;
            });
            if (res.point==p && p.isValid){
                ownerCache.pushFront(res);
                return res.owner;
            }
        }
        res.owner=pointOwner(SKeyVal.Any,p);
        res.point=p;
        if (p.isValid && input.ownerCacheSize>0){
            synchronized(this){
                ownerCache.filterInPlace(delegate bool(PointAndOwner po){
                    if (po.point==p){
                        res=po;
                        return false;
                    }
                    return true;
                });
                ownerCache.pushFront(res);
                if (ownerCache.length>input.ownerCacheSize) ownerCache.popEnd;
            }
        }
        return res.owner;
    }
    
    /// adds an extra observer that will be notified about the network changes
    void addObserver(ExplorationObserverI!(T) o){
        throw new Exploration("addition of observers not supported",__FILE__,__LINE__);
    }
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o){
        throw new Exploration("removal of observers not supported",__FILE__,__LINE__);
    }
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o){
        throw new Exploration("addition of explorers not supported",__FILE__,__LINE__);
    }
    /// removes the given explorer
    void rmExplorer(ExplorerI!(T)o){
        throw new Exploration("removal of explorers not supported",__FILE__,__LINE__);
    }
    
    /// reference position, this should be used just to create other ParticleSystem, never directly
    ParticleSys!(T) refPos(){
        return _refPos;
    }
    /// constraints for this system
    ConstraintI!(T) constraints(){
        return _constraints;
    }
    /// if the gradient is cheap to compute
    bool cheapGrad(){
        return properties["cheapGrad"]!=0;
    }
    /// local notification center
    NotificationCenter nCenter(){
        return _nCenter;
    }
    
    // values to define the point net topology
    mixin(extractFromProperties(propertiesList));
    
    /// creates a new point located at newPos in this silos, the point is not yet broadcasted
    /// not all silos might support creation of local points, use nextFreeSilos to get a silos
    /// thas supports it, use executeLocally to create, setup & publish a point...
    MainPointI!(T) newPointAt(DynPVector!(T,XType) newPos){
        throw new Exception("point creation not supported",__FILE__,__LINE__);
    }
    /// local point mainpoint (the real reference point)
    MainPointI!(T) mainPointL(Point){
        throw new Exception("no mainPoints owned",__FILE__,__LINE__);
    }
    /// publishes the local point given.
    /// Returns the pubblished point which might be different fron the argument if a collision did happen.
    MainPointI!(T) bcastPoint(MainPointI!(T)){
        throw new Exception("no mainPoints owned",__FILE__,__LINE__);
    }
    /// a local point (possibly a copy), is retained, and needs to be released (thus the create in the name)
    /// the time t is used to decide if a cached version can be used
    MainPointI!(T)createLocalPoint(Point p,Time t){
        if (hasKey(ownerOfPoint(p))){
            auto mp=localPoints[p];
            mp.retain;
            return mp;
        }
        CachedPoint!(T) cachedP;
        CachedPoint!(T) * pC;
        synchronized(localCache){
            pC= p in localCache;
            if (pC!is null){
                cachedP=*pC;
            }
        }
        if (pC is null){
            cachedP.mainPoint=pPool.getObj();
            cachedP.mainPoint.pos.checkX();
            cachedP.mainPoint._point=p;
            cachedP.mainPoint._gFlags|=GFlags.LocalCopy;
        } else if (cachedP.lastSync>t){
            cachedP.mainPoint.retain;
            return cachedP.mainPoint;
        }
        cachedP.lastSync=Clock.now;
        cachedP.mainPoint[]=mainPoint(ownerOfPoint(p),p);
        synchronized(localCache){
            localCache[p]=cachedP;
        }
        return cachedP.mainPoint;
    }
    /// drops a cached point (the point is not in use anymore)
    void dropCachedPoint(MainPointI!(T)p){
        synchronized(localCache){
            localCache.removeKey(p.point);
        }
    }

    /// an url that can be used to contact the silos core
    char[] silosCoreUrl(){
        return input.connectionUrl;
    }
}
