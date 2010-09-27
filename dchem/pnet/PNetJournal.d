/// keeps a journal of the changes
module dchem.pnet.PNetJournal;
import dchem.Common;
import dchem.pnet.PNetModels;

// object that keeps the journal of the computations done
class PNetJournalGen:ExplorationObserverGen{
    char[] journalDir;
    char[] journalFormat;
    bool logOtherEnergies;
    bool logOtherStart;
    bool logOtherPos;

    this(){}
    ExplorationObserverI!(Real) realObserver(LocalSilosI!(Real) silos){
        return new PNetJournal!(Real)(this,silos);
    }
    ExplorationObserverI!(LowP) lowPObserver(LocalSilosI!(LowP) silos){
        return new PNetJournal!(LowP)(this,silos);
    }
    mixin(serializeSome("dchem.PNetJournal",`
    journalDir: directory where to write the journal`))
}

// object that keeps the journal of the computations done
class PNetJournal(T):ExplorationObserverI{
    PNetJournalGen context;
    PNetSilos!(T) silos;
    char[] journalPath;
    Serializer jSerial;
    const char[] jVersion="PNetJournal v1.0";
    RLock serialLock;
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    HashSet!(Point) pointsToRm;
    
    struct JournalEntry{
        enum Kind:int{
            PointGrad, StartPoint,FinishedPoint,publishCollision,EvalE
        }
        Kind kind;
        uint flags;
        Point point;
        PSysWriter!(T) pos;
        SKey sKey;
        Real energy;
        mixin("dchem.JournalEntry!("~T.stringof~")","kind|flags|point|pos|energy");
        static JournalEntry PointGrad(Point p,PSysWriter!(T) pos){
            JournalEntry res;
            res.kind=Kind.PointGrad;
            res.flags=flags;
            res.point=point;
            return res;
        }
        static JournalEntry StartPoint(SKey sKey,Point p,PSysWriter!(T) pos,uint flags){
            JournalEntry res;
            res.kind=Kind.PointGrad;
            res.sKey=sKey;
            res.flags=flags;
            res.pos=pos;
            res.point=point;
            return res;
        }
        static JournalEntry FinishedPoint(Point p){
            JournalEntry res;
            res.kind=Kind.FinishedPoint;
            res.point=point;
            return res;
        }
        static JournalEntry publishCollision(Point p){
            JournalEntry res;
            res.kind=Kind.publishCollision;
            res.point=point;
            return res;
        }
        static JournalEntry EvalE(Point p,Real energy){
            JournalEntry res;
            res.kind=Kind.EvalE;
            res.point=point;
            res.energy=energy;
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
        }
        ClassMetaInfo getSerializationMetaInfo(){
            return metaI;
        }
        void serial(S s){
            bool skip=s.canDropDefaults();
            s.field(metaI[0],kind);
            if (!skip || kind==Kind.StartPoint) s.field(metaI[1],flags);
            s.field(metaI[2],point);
            if (!skip || pos !is null) s.field(metaI[3],pos);
            if (!skip || kind==Kind.EvalE) s.field(metaI[4],energy);
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
        this.context=c;
        this.silos=silos;
        this.serialLock=new RLock();
        long i=0;
        char[128] buf;
        auto arr=lGrowableArray(buf,0);
        OutputStream stream;
        while (i<100){
            const File.Style WriteUnique = {File.Access.Read, File.Open.New, File.Share.Read};
            arr(context.journalDir);
            writeOut(&arr.appendArr,i);
            try{
                stream=new File(arr.data,WriteUnique);
                this.journalPath=arr.takeData();
                break;
            } catch (Exception e){
                arr.clearData();
            } // should ignore only already present files... improve?
        }
        if (context.journalFormat=="json"){
            auto charSink=strDumper(stream);
            jSerial=new JsonSerializer!(T)(charSink);
        } else if (context.journalFormat=="sbin"){
            auto binSink=binaryDumper(stream);
            jSerial=new SbinSerializer(binSink);
        }
    }
    
    // ExplorationObserverI(T)
    
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(LocalSilosI!(T)silos,Point p,Real energy){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.EvalE(p,energy));
    }
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    /// sets inProgress to false
    void addGradEvalLocal(LocalSilosI!(T)silos,Point p,PSysWriter!(T) pSys){
        this.serialLock.lock();
        scope(exit){ this.serialLock.unlock(); }
        jSerial(JournalEntry.PointGrad(p,pSys));
    }
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void publishPoint(SKey owner,Point point,PSysWriter!(T) pos,T pSize,uint flags){
        if (context.logOtherStart || owner==silos.key){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.StartPoint(owner,point,
                ((context.logOtherPos || owner==silos.key)?pos:null),flags));
        }
    }
    /// finished exploring the given local point (remove it from the active points), bcasts finishedExploringPoint
    void finishedExploringPointLocal(Point p){
        if (!context.logOtherStart){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.FinishedPoint(p));
        }
    }
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(Point){
        if (context.logOtherStart){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.FinishedPoint(p));
        }
    }
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void publishCollision(SKey s,Point p){
        if (context.logOtherStart|| s==silos.key){
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            jSerial(JournalEntry.publishCollision(p));
        }
    }
}


