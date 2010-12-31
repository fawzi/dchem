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
import blip.io.StreamConverters;
import tango.io.device.File;
import tango.io.device.Conduit;
import dchem.input.RootInput;
import blip.io.EventWatcher:ev_time;

// object that keeps the journal of the computations done
class PNetJournalGen:ExplorationObserverGen{
    char[] journalDir;
    char[] journalFormat="sbin";
    bool logOtherEnergies;
    bool logOtherStart;
    bool logOtherPos;
    bool logNeighInfo;
    int maxJournalIds=10;

    this(){}
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos){
        return new PNetJournal!(Real)(this,silos);
    }
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos){
        return new PNetJournal!(LowP)(this,silos);
    }
    mixin(serializeSome("dchem.PNetJournal",`
    journalDir: directory where to write the journal
    journalFormat: to format of the journal: 'sbin' or 'json'
    logOtherEnergies: if the energies of non local points should be logged
    logOtherStart: if the start of non local points should be logged
    logOtherPos: if the positions of non local points should be logged
    logNeighInfo: if the updates of energies or gradients of the neighbors should be logged`));
    mixin myFieldMixin!();
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
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    HashSet!(Point) pointsToRm;
    
    static struct JournalEntry{
        enum Kind:int{
            PointGrad, StartPoint,FinishedPoint,PublishCollision,EvalE,NeighborHasEnergy,NeighborHasGradient
        }
        Kind kind;
        uint flags;
        Point point;
        PSysWriter!(T) pos;
        SKey sKey;
        Real energy;
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
            res.kind=Kind.PointGrad;
            res.sKey=sKey;
            res.flags=flags;
            res.pos=pos;
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
        static JournalEntry EvalE(Point p,Real energy){
            JournalEntry res;
            res.kind=Kind.EvalE;
            res.point=p;
            res.energy=energy;
            return res;
        }
        static JournalEntry NeighborHasEnergy(SKey s,Point p,Point[] neighbors,Real energy){
            JournalEntry res;
            res.kind=Kind.NeighborHasEnergy;
            res.point=p;
            res.energy=energy;
            res.neighs=neighbors;
            return res;
        }
        static JournalEntry NeighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, Real energy,LocalSilosI!(T) localSilos){
            JournalEntry res;
            res.kind=Kind.NeighborHasGradient;
            res.point=p.point;
            auto mp=p.mainPoint(localSilos);
            res.pos=pSysWriter(mp.pos);
            res.energy=energy;
            res.neighs=neighbors;
            return res;
        }
    
        static ClassMetaInfo metaI;
        static this(){
            metaI=ClassMetaInfo.createForType!(typeof(*this))("dchem.JournalEntry!("~T.stringof~")");
            metaI.addFieldOfType!(Kind)("kind","kind of the entry");
            metaI.addFieldOfType!(Point)("point","point");
            metaI.addFieldOfType!(PSysWriter!(T))("pos","position and forces");
            metaI.addFieldOfType!(Real)("energy","potential energy");
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
            if (!skip || kind==Kind.StartPoint) s.field(metaI[4],flags);
            if (!skip || kind==Kind.NeighborHasEnergy|| kind==Kind.NeighborHasEnergy) s.field(metaI[5],neighs);
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
        OutputStream stream;
        Exception lastE;
        for (int i=0;i<input.maxJournalIds;++i){
            const File.Style WriteUnique = {File.Access.Write, File.Open.New, File.Share.Read};
            arr(input.journalDir);
            if (input.journalDir.length>0 && input.journalDir[$-1]!='/') arr("/");
            arr(input.myFieldName);
            arr("-");
            writeOut(&arr.appendArr,i);
            arr(".log");
            try{
                stream=new File(arr.data,WriteUnique);
                this.journalPath=arr.takeData();
                lastE=null;
                break;
            } catch (Exception e){
                arr.clearData();
                lastE=e;
            } // should ignore only already present files... improve?
        }
        if (lastE!is null){
            throw new Exception("exception trying to open journal file",__FILE__,__LINE__,lastE);
        }
        if (input.journalFormat=="json"){
            auto charSink=strDumper(stream);
            jSerial=new JsonSerializer!(char)(charSink);
        } else if (input.journalFormat=="sbin"){
            auto binSink=binaryDumper(stream);
            jSerial=new SBinSerializer(binSink);
        }
        jSerial("journal start @ ");
        jSerial(ev_time());
    }
    
    // ExplorationObserverI(T)
    /// informs that a silos did a shutdown with the given speed (0, waits for pending points and calls finishers
    /// less than 0: does not stop, larger than 0 stops immediately (no finishers))
    void increaseRunLevel(SKey s,RunLevel rLevel){}
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey s,Point p,Real energy){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.EvalE(p,energy));
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(SKey s,Point p,PSysWriter!(T) pSys){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.PointGrad(p,pSys));
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
        }
    }
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(SKey s,Point p,SKey owner){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.FinishedPoint(p,owner));
    }
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(SKey s,Point p,SKey owner){
        if (input.logOtherStart){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.FinishedPoint(p,owner));
        }
    }
    /// informs silos s that source has done the initial processing of point p0,
    /// and p0 is now known and collision free
    void didLocalPublish(SKey s,Point p0,SKey source){
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey s,Point p){
        if (input.logOtherStart|| silos.hasKey(s)){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.PublishCollision(p));
        }
    }
    
    /// a neighbor point has calculated its energy (and not the gradient)
    /// neighbors should be restricted to s
    void neighborHasEnergy(SKey s,Point p,Point[] neighbors,Real energy){
        if (input.logNeighInfo){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.NeighborHasEnergy(s,p,neighbors,energy));
        }
    }
    /// the neighbor point p has calculated its gradient (and energy)
    /// neighbors should be restricted to silos
    void neighborHasGradient(SKey s,LazyMPLoader!(T)p, Point[] neighbors, Real energy){
        if (input.logNeighInfo){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.NeighborHasGradient(s,p,neighbors,energy,silos));
        }
    }
}


