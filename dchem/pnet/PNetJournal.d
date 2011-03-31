/// keeps a journal of the changes
module dchem.pnet.PNetJournal;
import dchem.Common;
import dchem.pnet.PNetModels;
import blip.serialization.Serialization;
import blip.parallel.smp.Wait;
import blip.container.MinHeap;
import dchem.pnet.MainPoint;
import blip.container.HashSet;
import dchem.input.WriteOut;
import blip.container.GrowableArray;
import blip.io.BasicIO;
import blip.io.FileStream;
import dchem.input.RootInput;
import blip.io.EventWatcher:ev_time,ev_tstamp;
import blip.parallel.mpi.Mpi;

// object that keeps the journal of the computations done
class PNetJournalGen:ExplorationObserverGen{
    char[] fileBasePath="journal";
    char[] journalFormat="sbin";
    bool logOtherEnergies;
    bool logOtherStart;
    bool logOtherPos;
    bool logNeighInfo;
    bool flushEachEntry=false;

    this(){}
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos){
        return new PNetJournal!(Real)(this,silos);
    }
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos){
        return new PNetJournal!(LowP)(this,silos);
    }
    mixin(serializeSome("PNetJournal",`object that keeps a journal of the important changes (migth be used to restart an exploration)`,
    `fileBasePath: start of the journal filename (defaults to journal)
    journalFormat: to format of the journal: 'sbin' or 'json'
    logOtherEnergies: if the energies of non local points should be logged
    logOtherStart: if the start of non local points should be logged
    logOtherPos: if the positions of non local points should be logged
    logNeighInfo: if the updates of energies or gradients of the neighbors should be logged
    flushEachEntry: if each entry should be immediately flushed (false)`));
    mixin myFieldMixin!();
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        if (journalFormat!="sbin" && journalFormat!="json"){
            res=false;
            dumper(s)("invalid journalFormat '")(journalFormat)("', expected json or sbin in field '")(myFieldName)("'\n");
        }
        return res;
    }
}

// object that keeps the journal of the computations done
class PNetJournal(T):ExplorationObserverI!(T){
    PNetJournalGen input;
    LocalSilosI!(T) silos;
    char[] journalPath;
    Serializer jSerial;
    const char[] jVersion="PNetJournal v1.0";
    RLock serialLock;
    string name(){
        return "PNetJournal_"~input.myFieldName;
    }
    
    static struct JournalEntry{
        enum Kind:int{
            None,PointGrad, StartPoint,FinishedPoint,PublishCollision,EvalE,NeighborHasEnergy,
            NeighborHasGradient,PublishedLocalPoint
        }
        Kind kind;
        uint flags;
        Point point;
        PSysWriter!(T) pos;
        SKey sKey;
        Real energy;
        Real energyError;
        Point[] neighs;

        static JournalEntry PointGrad(Point p,PSysWriter!(T) pos){
            JournalEntry res;
            res.kind=Kind.PointGrad;
            res.point=p;
            res.pos=pos;
            return res;
        }
        static JournalEntry StartPoint(SKey sKey,Point p,PSysWriter!(T) pos,uint flags){
            JournalEntry res;
            res.kind=Kind.StartPoint;
            res.sKey=sKey;
            res.flags=flags;
            res.pos=pos;
            res.point=p;
            return res;
        }
        static JournalEntry PublishedLocalPoint(SKey sKey,Point p){
            JournalEntry res;
            res.kind=Kind.PublishedLocalPoint;
            res.sKey=sKey;
            res.point=p;
            return res;
        }
        static JournalEntry FinishedPoint(Point p,SKey owner){
            JournalEntry res;
            res.kind=Kind.FinishedPoint;
            res.sKey=owner;
            res.point=p;
            return res;
        }
        static JournalEntry PublishCollision(Point p){
            JournalEntry res;
            res.kind=Kind.PublishCollision;
            res.point=p;
            return res;
        }
        static JournalEntry EvalE(Point p,Real energy,Real energyError){
            JournalEntry res;
            res.kind=Kind.EvalE;
            res.point=p;
            res.energy=energy;
            res.energyError=energyError;
            return res;
        }
        static JournalEntry NeighborHasEnergy(SKey s,Point[] neighbors,PointEMin energy){
            JournalEntry res;
            res.kind=Kind.NeighborHasEnergy;
            res.point=energy.point;
            res.energy=energy.energy;
            res.neighs=neighbors~energy.minimum;
            res.flags=energy.id;
            return res;
        }
        static JournalEntry NeighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin energy,LocalSilosI!(T) localSilos){
            JournalEntry res;
            res.kind=Kind.NeighborHasGradient;
            res.point=p.point;
            auto mp=p.mainPoint(localSilos);
            res.pos=pSysWriter(mp.pos);
            res.energy=energy.energy;
            res.neighs=neighbors~energy.minimum;
            res.flags=energy.id;
            return res;
        }
    
        static ClassMetaInfo metaI;
        static this(){
            metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.JournalEntry!("~T.stringof~")",
                "encodes a change in the point network");
            metaI.addFieldOfType!(Kind)("kind","kind of the entry");
            metaI.addFieldOfType!(Point)("point","point");
            metaI.addFieldOfType!(PSysWriter!(T))("pos","position and forces");
            metaI.addFieldOfType!(Real)("energy","potential energy");
            metaI.addFieldOfType!(Real)("energyError","potential energy Error");
            metaI.addFieldOfType!(uint)("flags","flags of the point");
            metaI.addFieldOfType!(Point[])("neighs","neighbors of the point");
        }
        ClassMetaInfo getSerializationMetaInfo(){
            return metaI;
        }
        void serial(S)(S s){
            bool skip=s.canDropDefaults();
            s.field(metaI[0],kind);
            s.field(metaI[1],point);
            if ((!skip) || (!pos.isDummy())) s.field(metaI[2],pos);
            if (!skip || kind==Kind.EvalE || kind==Kind.NeighborHasGradient || kind==Kind.NeighborHasEnergy) s.field(metaI[3],energy);
            if (!skip || kind==Kind.EvalE) s.field(metaI[4],energyError);
            if (!skip || kind==Kind.StartPoint) s.field(metaI[5],flags);
            if (!skip || kind==Kind.NeighborHasEnergy|| kind==Kind.NeighborHasEnergy) s.field(metaI[6],neighs);
        }
        void serialize(Serializer s){
            serial(s);
        }
        void unserialize(Unserializer s){
            serial(s);
        }
    }
    
    /// creates a new journal logger
    this(PNetJournalGen c,LocalSilosI!(T) silos){
        this.input=c;
        this.silos=silos;
        this.serialLock=new RLock();
        char[128] buf;
        auto arr=lGrowableArray(buf,0);
        Exception lastE;
        if (lastE!is null){
            throw new Exception("exception trying to open journal file",__FILE__,__LINE__,lastE);
        }
        if (input.journalFormat=="json"){
            auto stream=silos.outfileForName(input.fileBasePath~"."~input.journalFormat~"Log",WriteMode.WriteAppend,StreamOptions.CharBase);
            jSerial=new JsonSerializer!(char)(stream);
        } else if (input.journalFormat=="sbin"){
            auto stream=silos.outfileForName(input.fileBasePath~"."~input.journalFormat~"Log",WriteMode.WriteAppend,StreamOptions.BinBase);
            jSerial=new SBinSerializer(stream);
        } else {
            assert(0);
        }
        jSerial(jVersion);
        jSerial("journal start @ ");
        jSerial(ev_time());
        jSerial(T.stringof);
        if (input.flushEachEntry) jSerial.flush();
    }
    
    // ExplorationObserverI(T)
    /// informs that a silos did a shutdown with the given speed (0, waits for pending points and calls finishers
    /// less than 0: does not stop, larger than 0 stops immediately (no finishers))
    void increaseRunLevel(SKey s,RunLevel rLevel){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial.flush();
        if (rLevel>RunLevel.Stopping) jSerial.close();
    }
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy,Real energyError){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.EvalE(p,energy,energyError));
        if (input.flushEachEntry) jSerial.flush();
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.PointGrad(p,pSys));
        if (input.flushEachEntry) jSerial.flush();
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void publishPoint(SKey s,SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
        PSysWriter!(T) dummy;
        if (input.logOtherStart || silos.hasKey(owner)){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.StartPoint(owner,point,
                ((input.logOtherPos || silos.hasKey(owner))?pos:dummy),flags));
            if (input.flushEachEntry) jSerial.flush();
        }
    }
    /// communicates that the given local point has been published
    /// flags: communicate doubleEval?
    void publishedLocalPoint(SKey s,Point point){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.PublishedLocalPoint(s,point));
        if (input.flushEachEntry) jSerial.flush();
    }
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(SKey s,Point p,SKey owner){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.FinishedPoint(p,owner));
        if (input.flushEachEntry) jSerial.flush();
    }
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point p,SKey owner){
        if (input.logOtherStart){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.FinishedPoint(p,owner));
            if (input.flushEachEntry) jSerial.flush();
        }
    }
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source,int level){
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey s,Point p){
        if (input.logOtherStart|| silos.hasKey(s)){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.PublishCollision(p));
            if (input.flushEachEntry) jSerial.flush();
        }
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to s
    void neighborHasEnergy(SKey s,Point[] neighbors,PointEMin energy){
        if (input.logNeighInfo){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.NeighborHasEnergy(s,neighbors,energy));
            if (input.flushEachEntry) jSerial.flush();
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, PointEMin energy){
        if (input.logNeighInfo){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.NeighborHasGradient(s,p,neighbors,energy,silos));
            if (input.flushEachEntry) jSerial.flush();
        }
    }
    /// should speculatively calculate the gradient? PNetSilos version calls addEnergyEvalLocal
    bool speculativeGradientLocal(SKey s,Point p,Real energy,Real energyError){ return false; }
    /// checks it local point is somehow invalid and should better be skipped
    bool shouldFilterLocalPoint(SKey s,Point p){ return false; }
}

// object that loads some journals
class PNetJournalLoaderGen:SilosWorkerGen{
    char[][] filePaths=[];
    char[] journalFormat="sbin";
    char[] journalAccuracy="";
    bool parallelLoad=true;

    this(){}
    SilosWorkerI!(Real) silosWorkerReal(){
        return new PNetJournalLoader!(Real)(this);
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        return new PNetJournalLoader!(LowP)(this);
    }
    mixin(serializeSome("PNetJournalLoader",`Loads one or more journals`,
    `filePaths: files to load
    journalFormat: to format of the journal: 'sbin' or 'json'
    journalAccuracy: Real or LowP, by default the current accuracy of the silos
    parallelLoad: if when there are several files they should be loaded in parallel (true)`));
    mixin myFieldMixin!();
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        if (journalFormat!="sbin" && journalFormat!="json"){
            res=false;
            dumper(s)("invalid journalFormat '")(journalFormat)("', expected json or sbin in field '")(myFieldName)("'\n");
        }
        return res;
    }
}

/// actual journal loader worker
class PNetJournalLoader(T):SilosWorkerI!(T){
    PNetJournalLoaderGen input;
    LocalSilosI!(T) silos;
    
    this(PNetJournalLoaderGen input){
        this.input=input;
    }
    void workOn(LocalSilosI!(T)silos){
        switch(input.journalAccuracy){
        case("Real"):
            loadWithPrecision!(Real)(1);
            break;
        case("LowP"):
            loadWithPrecision!(LowP)(1);
            break;
        case(""):
            loadWithPrecision!(T)(1);
            break;
        default:
            assert(0);
        }
    }

    void loadWithPrecision(V)(int tag){
        alias PNetJournal!(V).JournalEntry.Kind Kind;
        auto paraEnv=silos.paraEnv();
        Point[Point] remappedPoints;
        Unserializer[] jUnserials;
        PNetJournal!(V).JournalEntry[] entries;
        char[][] paths;
        { // start all files
            size_t i=paraEnv.myRank;
            if ((!input.parallelLoad)&&paraEnv.myRank!=0) i=input.filePaths.length;
            while(i<input.filePaths.length){
                try{
                    auto f=infile(input.filePaths[i]);
                    Unserializer jUnserial;
                    if (input.journalFormat=="json"){
                        auto reader=infileStr(input.filePaths[i]);
                        jUnserial=new JsonUnserializer!(char)(reader);
                    } else if (input.journalFormat=="sbin"){
                        auto reader=infileBin(input.filePaths[i]);
                        jUnserial=new SBinUnserializer(reader);
                    } else {
                        assert(0);
                    }
                    char[] vers;
                    jUnserial(vers);
                    if (vers!=PNetJournal!(V).jVersion){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("different version string for journal loaded from file ")(input.filePaths[i])(":")(vers);
                        }),__FILE__,__LINE__);
                    }
                    char[] startStr;
                    jUnserial(startStr);
                    ev_tstamp tm;
                    jUnserial(tm);
                    char[] accuracy;
                    jUnserial(accuracy);
                    if (accuracy!=V.stringof){
                        throw new Exception("unexpected accuracy:"~accuracy~" vs "~V.stringof,__FILE__,__LINE__);
                    }
                    PNetJournal!(V).JournalEntry jEntry;
                    jUnserial(jEntry);
                    jUnserials~=jUnserial;
                    entries~=jEntry;
                    paths~=input.filePaths[i];
                } catch (Exception e){
                    silos.logMsg(delegate void(CharSink s){
                        dumper(s)("exception loading file ")(input.filePaths[i])(", trying to ignore:")(e);
                    });
                }
                if (input.parallelLoad) {
                    i+=paraEnv.dim();
                } else {
                    i+=1;
                }
            }
        }
pointLoop:while(1){
            Point minPoint=Point(0);
            size_t idxMinPoint=jUnserials.length;
            for (size_t i=0;i<jUnserials.length;++i){ // jUnserials.length might change during loop!!!
noPointStart:   while(1){ // proceed until next point
                    auto entry=&(entries[i]);
                    Point localP=entry.point;
                    auto lp=localP in remappedPoints;
                    if (lp!is null){
                        localP=*lp;
                    }
                    switch(entry.kind){
                        case Kind.StartPoint:
                            if (minPoint.data==0 || minPoint.data>entry.point.data){
                                minPoint=entry.point;
                                idxMinPoint=i;
                            }
                            break noPointStart;
                        case Kind.PointGrad:
                            static if (is(V==T)){
                                silos.addGradEvalLocal(silos.ownerOfPoint(localP),localP,entry.pos);
                            } else {
                                assert(0,"not implemented");
                            }
                            break;
                        case Kind.PublishCollision:
                            /// no point with StartPoint should have a collision, so we should be able to ignore it
                            break;
                        case Kind.EvalE:
                            silos.addEnergyEvalLocal(silos.ownerOfPoint(localP),localP,entry.energy,entry.energyError);
                            break;
                        default:
                    }
                    try{
                        jUnserials[i](*entry);
                    } catch (Exception e){
                        silos.logMsg(delegate void(CharSink s){
                            dumper(s)("stopping loading file ")(input.filePaths[i])(" due to:")(e);
                        });
                        entries[i]=entries[$-1];
                        jUnserials[i]=jUnserials[$-1];
                        paths[i]=paths[$-1];
                        entries=entries[0..$-1];
                        jUnserials=jUnserials[0..$-1];
                        paths=paths[0..$-1];
                        if (i==jUnserials.length) break noPointStart;
                    }
                }
            }
            auto points=new Point[](paraEnv.dim);
            points[paraEnv.myRank]=minPoint;
            if (input.parallelLoad){
                mpiAllGatherT(paraEnv,minPoint,points,tag);
            }
            Point pAtt=Point(0);
            size_t rankAtt=paraEnv.dim();
            foreach (iDim,p;points){
                if (p.data<pAtt.data||(pAtt.data==0&&p.data!=0)){
                    pAtt=p;
                    rankAtt=iDim;
                }
            }
            if (rankAtt==-1) break pointLoop;
            auto nSilos=silos.nextSilosRank(tag);
            PSysWriter!(V) pos;
            if (rankAtt!=nSilos){
                if (paraEnv.myRank==rankAtt){
                    auto sender=paraEnv[nSilos].sendTag(tag);
                    sender(entries[idxMinPoint].pos);
                    sender.close();
                }
                if (minPoint!=pAtt){
                    auto rcvr=paraEnv[rankAtt].recvTag(tag);
                    rcvr(pos);
                    rcvr.close();
                }
            } else {
                pos=entries[idxMinPoint].pos;
            }
            if (rankAtt==nSilos){
                auto posN=silos.refPos.dup();
                static if (is(V==T)){
                    posN[]=pos;
                } else {
                    assert(0,"not implemented");
                }
                auto newP=silos.newPointAt(posN.dynVars.x,pAtt);
                if (!isNullT(posN.dynVars.dx)){
                    silos.addGradEvalLocal(silos.ownerOfPoint(newP.point),newP.point,pSysWriter(posN));
                }
                minPoint=newP.point;
                posN.release();
            }
            mpiBcastT(paraEnv,minPoint,nSilos,tag);
            if (minPoint!=pAtt){
                remappedPoints[pAtt]=minPoint;
            }
            if (rankAtt==paraEnv.myRank){
                entries[idxMinPoint].kind=Kind.None;
            }
        }
    }
}
