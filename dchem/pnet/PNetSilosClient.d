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
import dchem.pnet.MainPoint;
import blip.io.EventWatcher;
import blip.container.HashMap;
import blip.container.GrowableArray;
import blip.parallel.mpi.MpiModels;

/// help structure 
struct CachedPoint(T){
    MainPoint!(T) mainPoint;
    ev_tstamp lastSync;
}

/// parameters for PNetSilosClient
class PNetSilosClientGen:InputElement {
    char[] connectionUrl;
    long ownerCacheSize=10;
    mixin(serializeSome("dchem.PNetSilosClient",`
    connectionUrl: the url to use to get the connection to the silos
    ownerCacheSize: amount of points whose owner is cached (10)`));
    mixin printOut!();
    mixin myFieldMixin!();
    this(){ }
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

class PNetSilosClient(T): LocalSilosI!(T){
    PNetSilosClientGen input;
    PNetSilosI!(T) connection;
    Real[char[]] properties;
    RandomSync _rand;
    CharSink log;
    Deque!(PointAndOwner) ownerCache;
    ParticleSys!(T) _refPos;
    ConstraintI!(T) _constraints;
    HashMap!(Point,CachedPoint!(T)) localCache;
    NotificationCenter _nCenter;
    CalculationContext ctx;
    PoolI!(MainPoint!(T)) pPool;
    
    MainPoint!(T)allocPoint(PoolI!(MainPoint!(T))p){
        return new MainPoint!(T)(this,p);
    }
    
    this(PNetSilosClientGen input, CalculationContext ctx, PNetSilosI!(T) connection=null,NotificationCenter nCenter=null,PoolI!(MainPoint!(T))pPool=null,CharSink log=null){
        this.log=log;
        if (log is null) this.log=sout.call;
        this.input=input;
        this.ctx=ctx;
        this.pPool=pPool;
        this.ownerCache=new Deque!(PointAndOwner)();
        this.localCache=new HashMap!(Point,CachedPoint!(T))();
        if (pPool is null) this.pPool=cachedPool(&this.allocPoint);
        if (ctx is null){
            assert(0,"CalculationContext should not be null");
            //input.evaluator.contentT!(Method)().getCalculator(true,[]);
        }
        if (_refPos is null){
            switch(ctx.activePrecision()){
            case Precision.Real:
                _refPos=ctx.refPSysReal.dupT!(T)(PSDupLevel.All);
                break;
            case Precision.LowP:
                _refPos=ctx.refPSysLowP.dupT!(T)(PSDupLevel.All);
                break;
            default:
                assert(0,"unknown activePrecision");
            }
        }
        ConstraintGen constraintGen=ctx.constraintGen;
        if (constraintGen !is null){
            _constraints=constraintT!(T)(constraintGen,_refPos);
        }
        if (_constraints is null) _constraints=new NoConstraint!(T)();
        
        this.connection=connection;
        if (connection is null){
            logMsg(delegate void(CharSink s){
                dumper(s)("connectionUrl:")(input.connectionUrl);
            });
            this.connection=ProtocolHandler.proxyForUrlT!(PNetSilosI!(T))(input.connectionUrl);
        }
        this._nCenter=nCenter;
        if (nCenter is null){
            this._nCenter=new NotificationCenter();
        }
        _rand=new RandomSync();
        properties=this.propertiesDict(SKeyVal.Any);
    }
    // ExplorationObserverI
    
    string name(){
        return "PNetSilosClient_"~input.myFieldName;
    }
    /// updates a pending operation status, should remove the operation when finished
    void updateEvalStatus(SKey owner,char[] opId, ModifyEvalOp!(T) op,EvalOp!(T).Status newStat){
        connection.updateEvalStatus(owner,opId,op,newStat);
    }
    /// this should be called by the master process sending it to SKey.All
    /// to receive a new operation to do
    void prepareNextOp(SKey s, int tag){
        connection.prepareNextOp(s,tag);
    }
    /// inserts an operation to execute into the server (target should be SKeyVal.Master)
    void addEvalOp(SKey s,EvalOp!(T)op,bool incrementNPending){
        connection.addEvalOp(s,op,incrementNPending);
    }
    /// stops a silos
    /// at the moment there is no support for dynamic adding/removal of silos, (worth adding???)
    /// so s should be only SKeyVal.All
    void increaseRunLevel(SKey s,RunLevel level){
        if (s!= SKeyVal.All) throw new Exception("should be called only with SKeyVal.All");
        connection.increaseRunLevel(s,level);
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
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin eAndMin){
        connection.neighborHasEnergy(s,neighbors,eAndMin);
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin eAndMin){
        connection.neighborHasGradient(s,p,neighbors,eAndMin);
    }
    
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point p,SKey owner){
        connection.finishedExploringPoint(s,p,owner);
    }
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source,int level){
        connection.didLocalPublish(s,p0,source,level);
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
    
    /// returns the next operation to execute
    EvalOp!(T) getNextOp(SKey s){
        version(TrackPNet) logMsg1("trying to getNextOp");
        auto res=connection.getNextOp(s);
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("nextOp is ")(res);
        });
        return res;
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
    
    /// energy for the local points and the minimum they reach (NAN if not yet known)
    PointEMin[] energyForPointsLocal(SKey s,Point[] pts,PointEMin[] ens){
        return connection.energyForPointsLocal(s,pts,ens);
    }
    /// energy for the points and the minimum they reach (NAN if not yet known)
    PointEMin[] energyForPoints(SKey s,Point[] pts,PointEMin[] ens){
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
    int nextSilosRank(int tag){
        assert(0,"implemented only in core silos");
    }
    void activateExplorer(SKey key,char[] name){
        return connection.activateExplorer(key,name);
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
        sinkTogether(log,delegate void(CharSink s){
            dumper(s)("<PNetSilosClientLog field=\"")(input.myFieldName)("\" addr=")(cast(void*)this)
                (" time=")(ev_time())(" task=\"")(taskAtt.val)("\" >\n  ");
            indentWriter("  ",s,writer);
            dumper(s)("\n</PNetSilosClientLog>\n");
        });
    }
    /// writes out a log message
    void logMsg1(char[]msg){
        sinkTogether(log,delegate void(CharSink s){
            dumper(s)("<PNetSilosClientLog field=\"")(input.myFieldName)("\" addr=")(cast(void*)this)
                (" time=")(ev_time())(" task=\"")(taskAtt.val)("\" >\n  ");
            sinkIndented("  ",s,msg);
            dumper(s)("\n</PNetSilosClientLog>\n");
        });
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
                if (ownerCache.length>input.ownerCacheSize) {
                    PointAndOwner dropped;
                    ownerCache.popBack(dropped);
                }
            }
        }
        return res.owner;
    }
    
    /// adds an extra observer that will be notified about the network changes
    void addObserver(ExplorationObserverI!(T) o){
        throw new Exception("addition of observers not supported",__FILE__,__LINE__);
    }
    /// removes the given observer if present
    void rmObserver(ExplorationObserverI!(T) o){
        throw new Exception("removal of observers not supported",__FILE__,__LINE__);
    }
    /// adds an extra explorer that can generate new points
    void addExplorer(ExplorerI!(T)o){
        throw new Exception("addition of explorers not supported",__FILE__,__LINE__);
    }
    /// removes the given explorer
    bool rmExplorerNamed(string){
        return false; // throw??
    }
    /// notify observers, the operation should not raise (or the whole program stops)
    void notifyLocalObservers(void delegate(ExplorationObserverI!(T))n){
        // no local observers
    }
    /// adds a component
    void addComponent(SilosComponentI!(T)component){
        throw new Exception("addition of components not supported",__FILE__,__LINE__); // change???
    }
    /// removes the component with the given name (stopping it)
    bool rmComponentNamed(string name){
        return false; // throw??
    }
    /// returns the current components
    HashMap!(string,SilosComponentI!(T)) components(){
        return null; // throw?
    }
    
    /// linear communicator (valid only inside the real silos, not in the clients)
    LinearComm paraEnv(){
        return null;
    }
    /// reference position, this should be used just to create other ParticleSystem, never directly
    ParticleSys!(T) refPos(){
        return _refPos;
    }
    /// constraints for this system
    ConstraintI!(T) constraints(){
        return _constraints;
    }
    /// how eagerly the gradient is calculated
    GradEagerness gradEagerness(){
        return cast(GradEagerness)cast(uint)properties["gradEagerness"];
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
    MainPointI!(T) newPointAt(DynPVector!(T,XType) newPos,Point proposedPoint){
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
    MainPointI!(T)createLocalPoint(Point p,ev_tstamp t){
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("creating new local point ")(p);
        });
        CachedPoint!(T) cachedP;
        CachedPoint!(T) * pC;
        synchronized(localCache){
            pC= p in localCache;
            if (pC!is null){
                cachedP=*pC;
            }
        }
        if (pC !is null && cachedP.lastSync>t){
            cachedP.mainPoint.retain;
            version(TrackPNet) logMsg(delegate void(CharSink s){
                dumper(s)("will return ")(cachedP.mainPoint);
            });
            return cachedP.mainPoint;
        }
        cachedP.mainPoint=pPool.getObj();
        cachedP.mainPoint.pos.checkX();
        cachedP.mainPoint._point=p;
        cachedP.lastSync=ev_time();
        cachedP.mainPoint[]=mainPoint(ownerOfPoint(p),p);
        cachedP.mainPoint._gFlags|=GFlags.LocalCopy;
        synchronized(localCache){
            localCache[p]=cachedP;
        }
        version(TrackPNet) logMsg(delegate void(CharSink s){
            dumper(s)("will return ")(cachedP.mainPoint);
        });
        cachedP.mainPoint.retain;
        return cachedP.mainPoint;
    }
    /// drops a cached point (the point is not in use anymore)
    void dropCachedPoint(MainPointI!(T)p){
        synchronized(localCache){
            auto lP=p.point in localCache;
            if (lP!is null && lP.mainPoint is p){
                p.release;
            }
            localCache.removeKey(p.point);
        }
    }
    /// registers a pending operation (sets id,...)
    void registerPendingOp(EvalOp!(T)op){
        assert(0,"no pending ops in silos client");
    }

    /// an url that can be used to contact the silos core
    char[] silosCoreUrl(){
        return input.connectionUrl;
    }
    
    int opApply(int delegate(ref Point,ref MainPointI!(T) el)loopBody){ return 0; }
}
