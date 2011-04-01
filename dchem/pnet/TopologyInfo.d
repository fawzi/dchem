/// keeps info on minima, boundary regions and attractors
module dchem.pnet.TopologyInfo;
import dchem.Common;
import dchem.input.RootInput;
import dchem.pnet.PNetModels;
import dchem.sys.ParticleSys;
import dchem.pnet.PointEvalOp;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.io.BasicStreams;
import dchem.pnet.DirArray;
import blip.io.FileStream;
import blip.container.Set;
import blip.container.GrowableArray;
import blip.util.NotificationCenter;
import blip.core.Variant;
import dchem.pnet.MainPoint;
import blip.core.stacktrace.StackTrace;
import dchem.calculator.FileCalculator;
import dchem.calculator.Calculator;
import blip.core.Array;
import blip.container.HashMap;
import blip.core.Traits;
import dchem.pnet.EmptyObserver;
import dchem.pnet.DenseLocalPointArray;
import blip.container.BatchedGrowableArray;
import blip.container.MinHeap;
import blip.container.AtomicSLink;
import blip.math.IEEE;
import blip.math.Math;
import blip.parallel.smp.Wait;

/// Information about a point
struct PointTopoInfo{
    Point minimum;
    ulong id;
    BoundaryPointInfo *border;
    mixin(serializeSome("PointTopoInfo","min/border information about a point","minimum|id"));
    mixin printOut!();
}

/// two minima (always ordered so min1<=min2)
struct Mins{
    Point min1;
    Point min2;
    mixin(serializeSome("","ordered couple of minima","min1|min2"));
    mixin printOut!();
    
    static Mins opCall(Point min1,Point min2){
        Mins res;
        if (min1.data<min2.data){
            res.min1=min1;
            res.min2=min2;
        } else {
            res.min1=min2;
            res.min2=min1;
        }
        return res;
    }
    equals_t opEquals(Mins m2){
        return min1==m2.min1 && min2==m2.min2;
    }
    hash_t toHash(){
        static assert(Mins.sizeof%size_t.sizeof==0);
        return rt_hash_block(cast(size_t*)cast(void*)this,Mins.sizeof/size_t.sizeof);
    }
    int opCmp(Mins m2){
        if (min1.data<m2.min1.data){
            return -1;
        } else if (min1.data==m2.min1.data){
            return ((min2.data<m2.min2.data)?-1:((min2.data==m2.min2.data)?0:1));
        } else {
            return 1;
        }
    }
}

/// two points at the border between two minima
struct BoundaryPointInfo{
    Point p1; // the smaller of the two points
    Point p2; // the larger of the two points
    int dir12;
    int dir21;
    Real cartesianDist12;
    Real cartesianDist21;
    Real energy;
    bool shortest1;
    bool shortest2;
    BoundaryPointInfo *next;
    BoundaryPointInfo *next2;
    
    Point otherP(Point p){
        if (p==p1){
            return p2;
        } else if (p==p2){
            return p1;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    BoundaryPointInfo* nextP(Point p){
        if (p==p1){
            return next;
        } else if (p==p2){
            return next2;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    
    int dirP(Point p){
        if (p==p1){
            return dir12;
        } else if (p==p2){
            return dir21;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    void dirPSet(int v,Point p){
        if (p==p1){
            dir12=v;
        } else if (p==p2){
            dir21=v;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    
    bool shortestP(Point p){
        if (p==p1){
            return shortest1;
        } else if (p==p2){
            return shortest2;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    void shortestPSet(bool v,Point p){
        if (p==p1){
            shortest1=v;
        } else if (p==p2){
            shortest2=v;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }

    Real cartesianDistP(Point p){
        if (p==p1){
            return cartesianDist12;
        } else if (p==p2){
            return cartesianDist21;
        } else {
            throw new Exception("no point p in this boundary",__FILE__,__LINE__);
        }
    }
    
    mixin(serializeSome("BoundaryPointInfo","information in boundary point",
        "p1|p2|dir12|dir21|cartesianDist12|cartesianDist21|energy|shortest1|shortest2"));
    mixin printOut!();
    
    /// clears this struct
    void clear(){
        *this=BoundaryPointInfo.init;
    }
}

/// information on a border between two minima
class BorderInfo(T){
    Mins mins;
    BatchedGrowableArray!(BoundaryPointInfo) pointInfo;
    BoundaryPointInfo *freeList;
    SpecialCmpHeap!(size_t) tPoint;
    TopologyInfo!(T) topoInfo;
    mixin(serializeSome("","information of a border between two minima","mins"));
    mixin printOut!();
    this(){
        this(null,Mins(Point(0),Point(0)));
    }
    this(TopologyInfo!(T) topoInfo,Mins mins){
        this.topoInfo=topoInfo;
        this.mins=mins;
        pointInfo=new BatchedGrowableArray!(BoundaryPointInfo)();
        tPoint.cmpEls=&comparePIdx;
    }
    
    bool comparePIdx(size_t a,size_t b){
        return pointInfo[a].energy <= pointInfo[b].energy;
    }
    
    BoundaryPointInfo *transitionPoint(){
        size_t i;
        if (tPoint.first(i))
            return pointInfo.ptrI(i);
        return null;
    }
    
    BoundaryPointInfo *add(BoundaryPointInfo newInfo){
        BoundaryPointInfo * res=popFrom(freeList);
        if (res is null){
            res=pointInfo.appendElT(newInfo);
        } else {
            *res=newInfo;
        }
        auto oldP=transitionPoint();
        tPoint.push(pointInfo.ptr2Idx(res));
        if (transitionPoint!is oldP){
            logTPointChange(oldP,transitionPoint);
        }
        return res;
    }
    /// removes the given bInfo (that must be owned by this BorderInfo)
    void remove(BoundaryPointInfo *bInfo){
        bool shouldRemoveSelf=false;
        auto idx=pointInfo.ptr2Idx(bInfo);
        bool tPointChange=idx==tPoint.peek;
        if (!tPoint.remove(idx)){
            throw new Exception("did not find BoundaryPointInfo index in tPoint heap",__FILE__,__LINE__);
        }
        if (tPoint.length==0){
            synchronized(topoInfo){
                if (tPoint.length==0){
                    logTPointChange(bInfo,null);
                    auto p=mins in topoInfo.borders;
                    if (*p == this)
                        topoInfo.borders.removeKey(mins);
                }
            }
        }
        if (tPointChange) logTPointChange(pointInfo.ptrI(idx),transitionPoint);
        if (bInfo.p1.isValid){
            auto p1Info=topoInfo.pointTopoInfo.ptrI(bInfo.p1);
            BoundaryPointInfo ** bAtt=&(p1Info.border);
            bool found=false;
            while((*bAtt)!is null){
                if ((*bAtt) is bInfo){
                    *bAtt=(**bAtt).next;
                    found=true;
                    break;
                }
                if ((**bAtt).p1==bInfo.p1)
                    bAtt=&((**bAtt).next);
                else
                    bAtt=&((**bAtt).next2);
            }
            assert(found);
        }
        if (bInfo.p2.isValid){
            auto p2Info=topoInfo.pointTopoInfo.ptrI(bInfo.p2);
            BoundaryPointInfo ** bAtt=&(p2Info.border);
            bool found=false;
            while((*bAtt)!is null){
                if ((*bAtt) is bInfo){
                    *bAtt=(**bAtt).next2;
                    found=true;
                    break;
                }
                if ((**bAtt).p1==bInfo.p1)
                    bAtt=&((**bAtt).next);
                else
                    bAtt=&((**bAtt).next2);
            }
            assert(found);
        }
        version(TopologyInfoNoClear){} else {
            bInfo.clear();
            insertAt(freeList,bInfo);
        }
    }
    void removeAll(){
        auto v=pointInfo.view;
        if (v.length>0){
            logTPointChange(transitionPoint,null);
        }
        foreach(ref p;v){
            this.remove(&p);
        }
        freeList=null;
        pointInfo.deallocData();
    }
    /// merges borderInfo (b2 will be to removed)
    void mergeWith(BorderInfo!(T) b2){
        auto oldP=transitionPoint;
        BatchedGrowableArray!(BoundaryPointInfo) shortPointInfo=b2.pointInfo;
        BatchedGrowableArray!(BoundaryPointInfo) longPointInfo=pointInfo;
        if (pointInfo.length<b2.pointInfo.length){
            longPointInfo=b2.pointInfo;
            shortPointInfo=pointInfo;
        }
        auto v=shortPointInfo.view;
        foreach(ref p;v){
            if (p.p1.isValid && p.p2.isValid){
                BoundaryPointInfo * newVal=longPointInfo.appendElT(p);
                topoInfo.updateBInfo(&p,newVal);
            }
        }
        pointInfo=longPointInfo;
        b2.pointInfo=shortPointInfo; // tryDeleteT?
        if (transitionPoint!is oldP){
            logTPointChange(oldP,transitionPoint);
        }
    }
    bool updateEnergy(BoundaryPointInfo *bp){
        if (topoInfo.calcEnergy(bp)){
            auto idx=pointInfo.ptr2Idx(bp);
            tPoint.remove(idx);
            tPoint.push(idx);
            return true;
        }
        return false;
    }
    void logTPointChange(BoundaryPointInfo *oldTp,BoundaryPointInfo *newTp){
        if (oldTp !is null){
            topoInfo.logTp("remTp",oldTp);
        }
        if (newTp !is null){
            topoInfo.logTp("addTp",newTp);
        }
    }
}

/// information on a minimum
class MinInfo{
    Point minimum;
    Set!(Point) attractor;
    this(){
        this(Point(0));
    }
    this(Point min){
        this.minimum=min;
        this.attractor=new Set!(Point);
    }
    mixin(serializeSome("","information on a minimum","minimum"));
    mixin printOut!();
}

// minUpdates newMin
// attractor updates (for boundary points)
// neigh updates (for boundary points)

/// a loader of points (what is created by the input)
class TopologyInfoGen:ExplorationObserverGen{
    bool flushEachLine=true;
    char[] minInfoFile="mins.log";
    char[] tPointInfoFile="tPoints.log";
    EvalLog[] minLogs=[{targetFile:"minima.xyz",format:"xyz"}];
    EvalLog[] tPointLogs=[{targetFile:"tPoints.xyz",format:"xyz"}];
    EvalLog[] pathLogs=[{targetFile:"path.xyz",format:"xyz"}];
    Real convergenceLimit=0.866;
    
    mixin(serializeSome("TopologyInfo",`Writes out the special points found so far`,
    `minInfoFile: base path used for the file where the minima are logged, if empty no log is written (defaults to mins.log)
    tPointInfoFile: base path used for the file where the transition points are logged, if emty no log is written (defaults to tPoints.log)
    flushEachLine: if each line should be flushed (true)
    convergenceLimit: minimum cos value of the residual of self pointing gradients before minima are merged (set it larger than 1 to disable merging) (0.866)
    minLogs: what is logged for points that are suspected to be minima (minima.xyz)
    tPointLogs: what is logged for points that are suspected to be transition points (tPointLogs.xyz)
    pathLogs: what is logged for paths (*path.xyz)`));
    mixin printOut!();
    mixin myFieldMixin!();
    
    this(){}
    
    bool verify(CharSink s){
        bool res=true;
        res=res && EvalLog.checkLogs(s,minLogs,"minLogs",myFieldName);
        res=res && EvalLog.checkLogs(s,minLogs,"tPointLogs",myFieldName);
        res=res && EvalLog.checkLogs(s,minLogs,"pathLogs",myFieldName);
        return res;
    }
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real)silos){
        auto res=new TopologyInfo!(Real)(this,silos);
        return res;
    }
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP)silos){
        auto res=new TopologyInfo!(LowP)(this,silos);
        return res;
    }
}

class TopologyInfo(T):EmptyObserver!(T){
    LocalSilosI!(T) silos;
    TopologyInfoGen input;
    Callback *myCallback;
    OutStreamI minLog;
    OutStreamI tPointLog;
    DenseLocalPointArray!(PointTopoInfo) pointTopoInfo;
    HashMap!(Mins,BorderInfo!(T)) borders;
    HashMap!(Point,MinInfo) minima;
    RLock minimaLock; /// lock for the input.minLogs
    RLock tPointsLock; /// lock for the input.tPointLogs
    
    this(TopologyInfoGen input,LocalSilosI!(T)silos){
        this.input=input;
        silos=silos;
        this.minimaLock=new RLock();
        this.tPointsLock=new RLock();
        pointTopoInfo=new DenseLocalPointArray!(PointTopoInfo)(silos.localPointsKeys());
        if (input.minInfoFile.length>0) minLog=silos.outfileForName(input.minInfoFile,WriteMode.WriteAppend,StreamOptions.CharBase|StreamOptions.Sync);
        if (input.tPointInfoFile.length>0) tPointLog=silos.outfileForName(input.tPointInfoFile,WriteMode.WriteAppend,StreamOptions.CharBase|StreamOptions.Sync);
        workOn(silos);
    }
    
    void workOn(LocalSilosI!(T)silos){
        this.silos=silos;
        if (input.minInfoFile.length>0){
            this.minLog=silos.outfileForName(input.minInfoFile,WriteMode.WriteAppend,
                StreamOptions.CharBase|StreamOptions.Sync);
        }
        if (input.tPointInfoFile.length>0){
            this.tPointLog=silos.outfileForName(input.tPointInfoFile,WriteMode.WriteAppend,
                StreamOptions.CharBase|StreamOptions.Sync);
        }
        myCallback=silos.nCenter.registerCallback("attractorChange",
            &attractorChanged,Callback.Flags.ReceiveAll);
        //silos.addComponent(this); // avoid??
    }
    
    char[] name(){
        return "TopologyInfo_"~input.myFieldName;
    }
    char[] kind(){
        return "TopologyInfo";
    }
    void logTp(string op, BoundaryPointInfo* tp){
        if (tp is null) return;
        if (tPointLog!is null){
            sinkTogether(tPointLog.charSink,delegate void(CharSink s){
                PointTopoInfo *p1Info=pointTopoInfo.ptrI(tp.p1);
                PointTopoInfo *p2Info=pointTopoInfo.ptrI(tp.p2);
                auto mins=Mins(p1Info.minimum,p2Info.minimum);
                bool first=mins.min1==p1Info.minimum;
                dumper(s)(op)("\t ")(mins.min1.data)("\t ")((first?tp.p1.data:tp.p2.data))("\t ")
                    ((first?tp.p2.data:tp.p1.data))("\t ")(mins.min2.data)("\n");
            });
            if (input.flushEachLine) tPointLog.flush;
        }
        if (op=="addTp" && input.tPointLogs.length>0){
            for (int i=0;i<2;++i){
                MainPointI!(T) mp=silos.mainPointL(tp.p1);
                if (i==1) mp=silos.mainPointL(tp.p2);
                char[128] buf;
                auto arr=lGrowableArray(buf,0,GASharing.Local);
                writeOut(&arr.appendArr,mp.point.data);
                auto extRef=arr.takeData;
                foreach(l;input.tPointLogs){
                    auto f=outfile(l.targetFile,WriteMode.WriteAppend,StreamOptions.BinBase);
                    scope(exit){ f.flush(); f.close(); }
                    WriteOut.writeConfig(f,mp.pos,l.format,extRef);
                }
            }
        }
    }
    /// merges two minima
    void mergeMin(Point oldMin,Point newMin){
        logMin("remMin",newMin);
        bool shouldLogNewMin=false;
        synchronized(this){
            assert(oldMin in minima,"mergeMin with non known oldMin");
            auto oldMinInfo=minima[oldMin];
            oldMinInfo.minimum=newMin;
            Mins[128] buf;
            BorderInfo!(T)[128] buf2;
            auto toRm=lGrowableArray(buf,0);
            auto toAdd=lGrowableArray(buf2,0);
            foreach (mins,b;borders){
                if (mins.min1==oldMin || mins.min2==oldMin){
                    toRm(mins);
                    if (mins.min1==newMin || mins.min2==newMin){
                        b.removeAll();
                    } else {
                        logTp("remTp",b.transitionPoint());
                        if (mins.min1==oldMin){
                            b.mins=Mins(newMin,mins.min2);
                        } else {
                            b.mins=Mins(newMin,mins.min1);
                        }
                        auto biAtt=b.mins in borders;
                        if (biAtt!is null){
                            biAtt.mergeWith(b);
                        } else {
                            toAdd(b);
                            if (tPointLog!is null){
                                BoundaryPointInfo *tp=b.transitionPoint;
                                PointTopoInfo *p1Info=pointTopoInfo.ptrI(tp.p1);
                                PointTopoInfo *p2Info=pointTopoInfo.ptrI(tp.p2);
                                if (p1Info.minimum==oldMin) {
                                    p1Info.minimum=newMin;
                                } else if (p2Info.minimum==oldMin){
                                    p2Info.minimum=newMin;
                                } else assert(0);
                                logTp("addTp",b.transitionPoint);
                            }
                        }
                    }
                }
            }
            foreach (b;toRm.data){
                borders.removeKey(b);
            }
            foreach (b;toAdd.data){
                borders[b.mins]=b;
            }
            toAdd.deallocData();
            toRm.deallocData();
            auto newMinInfoP=newMin in minima;
            if(newMinInfoP !is null){
                auto newMinInfo= *newMinInfoP;
                if (newMinInfo.attractor.size<oldMinInfo.attractor.size){
                    foreach(p;newMinInfo.attractor){
                        oldMinInfo.attractor.add(p);
                    }
                    tryDeleteT(newMinInfo.attractor);
                    newMinInfo.attractor=oldMinInfo.attractor;
                    oldMinInfo.attractor=null;
                } else {
                    foreach(p;oldMinInfo.attractor){
                        newMinInfo.attractor.add(p);
                        (*pointTopoInfo.ptrI(p)).minimum=newMin;
                    }
                    tryDeleteT(oldMinInfo.attractor);
                    oldMinInfo.attractor=null;
                }
                minima.removeKey(oldMin);
            } else {
                minima[newMin]=oldMinInfo;
                minima.removeKey(oldMin);
                shouldLogNewMin=true;
            }
        }
        auto newMinMP=silos.mainPointL(newMin);
        auto oldMinMP=silos.mainPointL(oldMin);
        Real newMinE;
        ulong newMinId;
        synchronized(newMinMP){
            newMinE=newMinMP.energy;
            newMinId=newMinMP.attractor.id;
            if (shouldLogNewMin){
                logMin("addMin",newMin);
            }
        }
        Attractor attractor;
        bool didChangeAttractor=false;
        synchronized(oldMinMP){
            if (oldMinMP.attractor.minimum!=newMin){
                assert(newMinE<=attractor.energyThroughPoint || isNaN(attractor.energyThroughPoint));
                attractor=oldMinMP.attractor;
                attractor.throughPoint=newMin;
                attractor.energyThroughPoint=newMinE;
                attractor.idThroughPoint=newMinId;
                attractor.minimum=newMin;
                didChangeAttractor=true;
            }
        }
        if (didChangeAttractor) {
            oldMinMP.attractor=attractor;
        }
        if (minLog!is null){
            sinkTogether(minLog.charSink,delegate void(CharSink s){
                dumper(s)("minMerge\t ")(oldMin.data)("\t ")(newMin.data)("\n");
            });
            if (input.flushEachLine) minLog.flush();
        }
    }
    
    void logMin(string op,Point min){
        if (minLog!is null){
            sinkTogether(minLog.charSink,delegate void(CharSink s){
                dumper(s)(op)("\t ")(min)("\n");
            });
            if (input.flushEachLine) minLog.flush();
        }
        if (op=="addMin" && input.minLogs.length>0){
            MainPointI!(T) mp=silos.mainPointL(min);
            char[128] buf;
            auto arr=lGrowableArray(buf,0,GASharing.Local);
            writeOut(&arr.appendArr,mp.point.data);
            auto extRef=arr.takeData;
            foreach(l;input.minLogs){
                auto f=outfile(l.targetFile,WriteMode.WriteAppend,StreamOptions.BinBase);
                scope(exit){ f.flush(); f.close(); }
                WriteOut.writeConfig(f,mp.pos,l.format,extRef);
            }
        }
    }
    
    /// replaces oldBI with newBI (when the address changed)
    void updateBInfo(BoundaryPointInfo *oldBI,BoundaryPointInfo *newBI){
        // update p1
        BoundaryPointInfo **prev=&(pointTopoInfo.ptrI(oldBI.p1).border);
        bool found=false;
        while (*prev!is null){
            if ((*prev) is oldBI){
                (*prev)=newBI;
                found=true;
                break;
            }
            if ((**prev).p1==oldBI.p1){
                prev=&((**prev).next);
            } else {
                prev=&((**prev).next2);
            }
        }
        assert(found);
        // update p2
        prev=&(pointTopoInfo.ptrI(oldBI.p2).border);
        found=false;
        while (*prev!is null){
            if ((*prev) is oldBI){
                (*prev)=newBI;
                found=true;
                break;
            }
            if ((**prev).p1==oldBI.p2){
                prev=&((**prev).next);
            } else {
                prev=&((**prev).next2);
            }
        }
        assert(found);
    }
    /// minimum of a point changed
    void pointAttractorChanged(PointEMin newMin){
        auto point=newMin.point;
        PointTopoInfo *oldMin=pointTopoInfo.ptrI(point);
        if (newMin.id<=oldMin.id) return;
        bool checkB=false;
        synchronized(this){
            PointTopoInfo oldInfo=*oldMin;
            PointTopoInfo newInfo;
            newInfo.minimum=newMin.minimum;
            newInfo.id=newMin.id;
            newInfo.border=oldInfo.border;
            if (newInfo.id<=oldInfo.id) return;
            *oldMin=newInfo;
            if (oldMin.minimum == newMin.minimum) return;
            if (point in minima){ // minimum moved
                mergeMin(point,newMin.minimum);
            } else {
                if (oldMin.minimum in minima) minima[oldMin.minimum].attractor.remove(point);
                if (newMin.minimum in minima){
                    minima[newMin.minimum].attractor.add(point);
                } else {
                    auto newMinInfo=new MinInfo(newMin.minimum);
                    newMinInfo.attractor.add(point);
                    minima[newMin.minimum]=newMinInfo;
                    logMin("addMin",newMin.point);
                }
                if (oldMin.border !is null){
                    fixBoundaryMinChange(oldMin,point,oldMin.minimum,newMin.minimum);
                }
                checkB=true;
            }
        }
        if (checkB) checkForBoundary(point);
    }
    
    void attractorChanged(cstring notificationName,Callback* callback,Variant oldF){
        PointEMin* newMin=oldF.get!(PointEMin*)();
        pointAttractorChanged(*newMin);
    }
    
    /// called after a minimum of a single point with boundary has changed
    void fixBoundaryMinChange(PointTopoInfo *pInfo,Point point, Point oldMin, Point newMin){
        BoundaryPointInfo *bAtt=pInfo.border;
        while (bAtt!is null){
            bool first=true;
            Point other=bAtt.p2;
            if (other==point){
                other=bAtt.p1;
                first=false;
            }
            auto otherPInfo=pointTopoInfo[other];
            auto otherMin=otherPInfo.minimum;
            auto oldBorder=borders[Mins(oldMin,otherMin)];
            if (otherMin!=newMin){
                auto biAtt=Mins(newMin,otherMin) in borders;
                if (biAtt!is null){
                    biAtt.add(*bAtt);
                } else {
                    auto newBorder=new BorderInfo!(T)(this,Mins(newMin,otherMin));
                    borders[newBorder.mins]=newBorder;
                    newBorder.add(*bAtt);
                }
            }
            auto oldVal=bAtt;
            if (first){
                bAtt.p1=Point(0);
                bAtt=bAtt.next;
            } else {
                bAtt.p2=Point(0);
                bAtt=bAtt.next2;
            }
            oldBorder.remove(oldVal);
        }
    }
    
    static struct BPoint{
        size_t idx;
        Real distance;
    }
    /// checks if p is at a boundary (for non local points one should first check the closest point in each dir,and then ftech them at once)
    void checkForBoundary(Point p){
        MainPointI!(T) mp=silos.mainPointL(p);
        auto myInfo=pointTopoInfo.ptrI(p);
        BoundaryPointInfo[128] buf;
        auto newBorder=lGrowableArray(buf,0);
        synchronized(mp){
            if (isNaN(mp.energy)) return;
            Point oldP=p;
            PointAndDir[] neighs=mp.neighbors().data();
            DirDistances!(T)[] neighDists=mp.neighDistances().data();
            BPoint[Point] min2Dist;
            foreach (i,pAtt;neighs){
                if (isNaN(neighDists[i].energy)) continue;
                if (pAtt.point!=oldP){
                    oldP=pAtt.point;
                    auto pAttInfo=pointTopoInfo.ptrI(pAtt.point);
                    if (myInfo.minimum!=pAttInfo.minimum){
                        auto d=pAttInfo.minimum in min2Dist;
                        if (d!is null){
                            if (d.distance>neighDists[i].cartesianDist){
                                d.idx=i;
                                d.distance=neighDists[i].cartesianDist;
                            }
                        } else {
                            min2Dist[pAttInfo.minimum]=BPoint(i,neighDists[i].cartesianDist);
                        }
                    }
                }
            }
            auto nDirs=mp.ndirs;
            foreach(m,bP;min2Dist){
                auto dir=neighs[bP.idx].dir;
                if ((dir==0 || dir==nDirs) && neighs.length>bP.idx+1 && neighs[bP.idx+1].point==neighs[bP.idx].point){
                    dir=neighs[bP.idx+1].dir;
                }
                bool tooFar=false;
                if (dir==invalidDir || dir==nDirs) tooFar=true;
                foreach(i,pAtt;neighs){
                    if (pAtt.dir==dir && bP.idx!=i && (!isNaN(neighDists[i].energy)) 
                        && neighDists[i].cartesianDist<bP.distance){
                        tooFar=true;
                    }
                }
                if (!tooFar){
                    BoundaryPointInfo bNew;
                    if (neighs[bP.idx].point<p){
                        bNew.p1=neighs[bP.idx].point;
                        bNew.p2=p;
                        bNew.dir21=dir;
                        bNew.cartesianDist21=bP.distance;
                        bNew.shortest2=true;
                    } else {
                        bNew.p1=p;
                        bNew.p2=neighs[bP.idx].point;
                        bNew.dir12=dir;
                        bNew.cartesianDist12=bP.distance;
                        bNew.shortest1=true;
                    }
                    newBorder(bNew);
                }
            }
        }
        // has found border
        if (newBorder.length){
            synchronized(this){
                foreach (nb;newBorder.data){
                    maybeAddBoundaryPoint(nb);
                }
            }
        }
    }
    /// if not yet added and still to be added, adds BoundaryPoint nb
    bool maybeAddBoundaryPoint(BoundaryPointInfo nb){
        synchronized(this){
            auto nDirs=silos.ndirs;
            auto p1Data=pointTopoInfo.ptrI(nb.p1);
            auto p2Data=pointTopoInfo.ptrI(nb.p2);
            bool known=false;
            auto bNow=p1Data.border;
            while (bNow!is null){
                if (bNow.p1==nb.p1 && bNow.p2==nb.p2){
                    known=true;
                    break;
                }
            }
            if (known) return false;
            auto p1Mp=silos.mainPointL(nb.p1);
            auto p2Mp=silos.mainPointL(nb.p1);
            synchronized(p1Mp){
                PointAndDir[] neighs=p1Mp.neighbors().data();
                DirDistances!(T)[] neighDists=p1Mp.neighDistances().data();
                nb.dir12=invalidDir;
                nb.cartesianDist12=Real.max;
                nb.shortest1=true;
                foreach(i,p;neighs){
                    if (p.point==nb.p2){
                        if (nb.dir12==0 || nb.dir12==nDirs || nb.dir12==invalidDir){
                            nb.dir12=p.dir;
                            nb.cartesianDist12=neighDists[i].cartesianDist;
                        }
                    }
                }
                foreach(i,p;neighs){
                    if (p.dir==nb.dir12 && p.point!=nb.p2){
                        if (neighDists[i].cartesianDist<nb.cartesianDist12){
                            nb.shortest1=false;
                        }
                    }
                }
            }
            synchronized(p2Mp){
                PointAndDir[] neighs=p2Mp.neighbors().data();
                DirDistances!(T)[] neighDists=p2Mp.neighDistances().data();
                nb.dir21=invalidDir;
                nb.cartesianDist21=Real.max;
                nb.shortest2=true;
                foreach(i,p;neighs){
                    if (p.point==nb.p1){
                        if (nb.dir21==0 || nb.dir21==nDirs || nb.dir21==invalidDir){
                            nb.dir21=p.dir;
                            nb.cartesianDist21=neighDists[i].cartesianDist;
                        }
                    }
                }
                foreach(i,p;neighs){
                    if (p.dir==nb.dir21 && p.point!=nb.p1){
                        if (neighDists[i].cartesianDist<nb.cartesianDist21){
                            nb.shortest2=false;
                        }
                    }
                }
            }
            if (nb.shortest1||nb.shortest2){
                nb.energy=Real.init;
                auto mins=Mins(p1Data.minimum,p2Data.minimum);
                if (mins.min1!=mins.min2){
                    calcEnergy(&nb);
                    auto bAtt=mins in borders;
                    BoundaryPointInfo * newB;
                    if (bAtt){
                        newB=bAtt.add(nb);
                    } else {
                        auto bInfo=new BorderInfo!(T)(this,mins);
                        borders[mins]=bInfo;
                        newB=bInfo.add(nb);
                    }
                    newB.next=p1Data.border;
                    newB.next2=p2Data.border;
                    p1Data.border=newB;
                    p2Data.border=newB;
                    return true;
                }
            }
        }
        return false;
    }
    /// calculates the energy of the transition point through the given border point
    bool calcEnergy(BoundaryPointInfo *bp){
        auto p1Mp=silos.mainPointL(bp.p1);
        auto p2Mp=silos.mainPointL(bp.p2);
        Real e1,e2,proj1,proj2;
        Real oldEnergy=bp.energy;
        bool foundConvergence=false;
        synchronized(p1Mp){
            synchronized(p2Mp){
                e1=p1Mp.pos.dynVars.potentialEnergy;
                e2=p2Mp.pos.dynVars.potentialEnergy;
                PointAndDir[] neighs1=p1Mp.neighbors().data();
                DirDistances!(T)[] neighDists1=p1Mp.neighDistances().data();
                foreach(i,p;neighs1){
                    if (p.point==bp.p1){
                        proj1=neighDists1[i].minDirProj;
                        break;
                    }
                }
                PointAndDir[] neighs2=p2Mp.neighbors().data();
                DirDistances!(T)[] neighDists2=p2Mp.neighDistances().data();
                foreach(i,p;neighs2){
                    if (p.point==bp.p2){
                        proj2=neighDists2[i].minDirProj;
                        break;
                    }
                }
                if (isNaN(e1)){
                    if (isNaN(e2)){
                        bp.energy=Real.max;
                    } else {
                        bp.energy=e2;
                    }
                } else if (isNaN(e2)) {
                    bp.energy=e1;
                } else {
                    if (!isNaN(proj1) && !isNaN(proj2)){
                        if (proj1<=0 && proj2<=0){
                            Real a=e1,g=e2,b=proj1,f=-proj2;
                            Real c=3*g-f-3*a-2*b;
                            Real d=f-2*g+2*a+b;
                            auto A=3*d;
                            auto B=2*c;
                            auto C=b;
                            auto delta=B*B-4*A*C;
                            bp.energy=max(e1,e2);
                            if (delta>0){
                                auto x=(-b+sqrt(delta))/(2*A);
                                if (x>=0 && x<=1){
                                    auto nEnergy=a+x*(b+x*(c+x*d));
                                    if (nEnergy>bp.energy){
                                        bp.energy=nEnergy;
                                    } else {
                                        assert(0,"tmp check to remove");
                                    }
                                } else {
                                    assert(0,"tmp check to remove");
                                }
                            } else {
                                assert(0,"tmp check to remove");
                            }
                        } else if (proj1>0 && proj2>0){
                            // check for convergence
                            auto diff=p2Mp.pos.dynVars.x.dup();
                            diff.opBypax(p1Mp.pos.dynVars.x,-1,1);
                            auto deriv1=p1Mp.pos.dynVars.dVarStruct.emptyDx();
                            deriv1[]=0;
                            p1Mp.pos.addToTSpace!(T)(diff,deriv1);
                            auto deriv2=p2Mp.pos.dynVars.dVarStruct.emptyDx();
                            deriv2[]=0;
                            diff*=cast(T)-1;
                            p2Mp.pos.addToTSpace!(T)(diff,deriv2);
                            diff.giveBack();
                            auto minDirProj1=p1Mp.pos.dotInTSpace(deriv1,p1Mp.minDir);
                            auto minDirProj2=p2Mp.pos.dotInTSpace(deriv2,p2Mp.minDir);
                            auto rest1=p1Mp.pos.dynVars.mddx.dup();
                            silos.constraints().applyDR(p1Mp.pos,rest1);
                            rest1.opBypax(deriv1,-minDirProj1/p1Mp.minDirScale,1);
                            deriv1.giveBack();
                            auto rest2=p2Mp.pos.dynVars.mddx.dup();
                            silos.constraints().applyDR(p2Mp.pos,rest2);
                            rest2.opBypax(deriv2,-minDirProj2/p2Mp.minDirScale,1);
                            deriv2.giveBack();
                            auto r1N=rest1.norm2;
                            auto r2N=rest2.norm2;
                            auto r1r2=rest1.opDot(rest2);
                            if (r1N*r2N==0){
                                auto cosR=r1r2/(r1N*r2N);
                                if (cosR>input.convergenceLimit){
                                    foundConvergence=true;
                                }
                            }
                            rest1.giveBack();
                            rest2.giveBack();
                            bp.energy=max(e1,e2);
                        } else {
                            bp.energy=max(e1,e2);
                        }
                    } else {
                        bp.energy=max(e1,e2);
                    }
                }
            }
        }
        if (foundConvergence){
            synchronized(this){
                auto min1=pointTopoInfo[bp.p1].minimum;
                auto min2=pointTopoInfo[bp.p2].minimum;
                if (min1!=min2){
                    silos.logMsg(delegate void(CharSink s){
                        dumper(s)("detected convergent minima, merging ")(min1.data)(" and ")(min2.data)("\n");
                    });
                    auto min1E=silos.mainPointL(min1).energy;
                    auto min2E=silos.mainPointL(min2).energy;
                    if (min2E<min1E){
                        mergeMin(min1,min2);
                    } else {
                        mergeMin(min2,min1);
                    }
                }
            }
            return false;
        }
        return oldEnergy!=bp.energy;
    }
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to silos
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin energy){
        auto pInfo=pointTopoInfo.ptrI(energy.point);
        if (pInfo.id>=energy.id) return;
        if (pInfo.border!is null) {
            synchronized(this){
                auto bAtt=pInfo.border;
                while(bAtt!is null){
                    auto biAtt=borders[Mins(pointTopoInfo[bAtt.p1].minimum,pointTopoInfo[bAtt.p2].minimum)];
                    biAtt.updateEnergy(bAtt);
                    bAtt=bAtt.nextP(energy.point);
                }
            }
        }
        if (pInfo.minimum!=energy.minimum){
            pointAttractorChanged(energy);
        } else {
            bool checkB=false;
            foreach(p;neighbors){
                auto pTopo=pointTopoInfo.ptrI(p);
                if (energy.minimum!=pTopo.minimum){
                    if (pTopo.border!is null){
                        auto mp=silos.mainPointL(p);
                        Real cDist=Real.max;
                        int cDir=int.max;
                        synchronized(mp){
                            PointAndDir[] neighs=mp.neighbors().data();
                            DirDistances!(T)[] neighDists=mp.neighDistances().data();
                            foreach(i,pAtt;neighs){
                                if (p==pAtt.point){
                                    cDist=neighDists[i].cartesianDist;
                                    if (cDir==0||cDir==int.max||cDir==invalidDir){
                                        cDir=pAtt.dir;
                                        break;
                                    }
                                }
                            }
                        }
                        BoundaryPointInfo **bAtt=&pTopo.border;
                        while((*bAtt)!is null){
                            if ((**bAtt).dirP(p)==cDir && (**bAtt).cartesianDistP(p)<cDist && 
                                (**bAtt).p1!=energy.point && (**bAtt).p2!=energy.point)
                            {
                                auto otherPoint=(**bAtt).otherP(p);
                                auto borderAtt=Mins(pTopo.minimum,energy.minimum) in borders;
                                borderAtt.remove(*bAtt);
                            } else {
                                if ((**bAtt).p1==energy.point){
                                    bAtt=&((**bAtt).next);
                                } else {
                                    bAtt=&((**bAtt).next2);
                                }
                            }
                        }
                    }
                    checkB=true;
                    break;
                }
            }
            if (checkB){
                checkForBoundary(energy.point);
            }
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin energy){
        neighborHasEnergy(s,neighbors,energy);
    }
    
    override void increaseRunLevel(SKey s,RunLevel level){
        if (level>RunLevel.Stopping) stop();
    }
    
    void printPath(string baseName,BoundaryPointInfo*bp,EvalLog[]logs){
        if (bp is null || logs.length==0) return;
        Point[128] buf;
        auto toMin1=lGrowableArray(buf,0);
        auto p1Info=pointTopoInfo.ptrI(bp.p1);
        auto p2Info=pointTopoInfo.ptrI(bp.p2);
        auto p1m=p1Info.minimum;
        auto mins=Mins(p1m,p2Info.minimum);
        Point p;
        if (p1m==mins.min1){
            p=bp.p1;
        } else {
            p=bp.p2;
        }
        toMin1(p);
        {
            auto attractor=silos.mainPointL(p).attractor;
            while(attractor.throughPoint!=p){
                p=attractor.throughPoint;
                toMin1(p);
                attractor=silos.mainPointL(p).attractor;
            }
        }
        char[256] buf2;
        auto fName=lGrowableArray(buf2,0);
        scope(exit) {
            fName.deallocData;
        }
        fName(baseName);
        fName("-");
        writeOut(&fName.appendArr,mins.min1.data);
        fName("_");
        writeOut(&fName.appendArr,mins.min2.data);
        size_t lenBase=fName.length;
        foreach(l;logs){
            fName.growTo(lenBase);
            fName(l.targetFile);
            auto f=silos.outfileForName(fName.data,WriteMode.WriteAppend,StreamOptions.BinBase);
            scope(exit){ f.flush(); f.close(); }
            char[128] buf3;
            auto extRef=lGrowableArray(buf3,0);
            foreach_reverse(pAtt;toMin1.data){
                extRef.clearData();
                writeOut(extRef,pAtt.data);
                auto mp=silos.mainPointL(pAtt);
                WriteOut.writeConfig(f,mp.pos,l.format,extRef.data);
            }
            if (p1m==mins.min1){
                p=bp.p2;
            } else {
                p=bp.p1;
            }
            auto attractor=silos.mainPointL(p).attractor;
            while(true){
                extRef.clearData();
                writeOut(extRef,p.data);
                auto mp=silos.mainPointL(p);
                WriteOut.writeConfig(f,mp.pos,l.format,extRef.data);
                if (attractor.throughPoint==p) break;
                p=attractor.throughPoint;
                attractor=silos.mainPointL(p).attractor;
            }
        }
    }
    
    void stop(){
        Callback *cb;
        synchronized(this){
            if (myCallback){
                cb=myCallback;
                myCallback=null;
            }
        }
        if (cb!is null){
            cb.flags&= ~Callback.Flags.Resubmit;
            silos.nCenter.unregisterReceiveAllCallback("localPointChangedGFlags",cb);
            //silos.rmComponentNamed(this.name);
        }
        if (minLog!is null){
            minLog.close();
            minLog=null;
        }
        if (tPointLog!is null){
            tPointLog.close();
            tPointLog=null;
        }
        if (input.pathLogs.length>0){
            foreach(m,b;borders){
                printPath("finalPath",b.transitionPoint,input.pathLogs);
            }
        }
    }
}
