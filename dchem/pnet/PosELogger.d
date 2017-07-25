/// logs position and energy of the points
module dchem.pnet.PosELogger;
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
import blip.io.EventWatcher:ev_time;
import dchem.pnet.EmptyObserver;

// object that keeps the journal of the computations done
class PosELoggerGen:ExplorationObserverGen{
    char[] logFile="posELog.log";
    bool flushEachLine=false;
    char[] posFormat="f8.4";
    char[] eFormat="g12.6";
    this(){}
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos){
        return new PosELogger!(Real)(this,silos);
    }
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos){
        return new PosELogger!(LowP)(this,silos);
    }
    mixin(serializeSome("PosELogger",`Logs position and energy`,
    `logFile: filename used for the logs
    posFormat: format string to format the positions
    eFormat: format string to format the energy
    flushEachLine: if each line should be flushed (defaults to 0)`));
    mixin myFieldMixin!();
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        {
            LowP l;
            Real r;
            Exception lowPE,realE;
            try{
                writeOut(delegate void(char[] s){},r,posFormat);
            } catch (Exception e){
                realE=e;
            }
            try{
                writeOut(delegate void(char[] s){},l,posFormat);
            } catch (Exception e){
                lowPE=e;
            }
            if (lowPE!is null && realE !is null){
                realE.next=lowPE;
                dumper(s)("invalid posFormat in field ")(myFieldName)(", both Real and LowP cannot use it. Real failed with:")(realE)(" LowP:")(lowPE)("\n");
                res=false;
            } else if (lowPE!is null){
                dumper(s)("warning partially valid posFormat in field ")(myFieldName)(", it will fail if used with a low precision point\n");
            } else if (realE!is null){
                dumper(s)("warning partially valid posFormat in field ")(myFieldName)(", it will fail if used with a real precision point\n");
            }
        }
        {
            Real r;
            Exception realE;
            try{
                writeOut(delegate void(char[] s){},r,eFormat);
            } catch (Exception e){
                dumper(s)("invalid eFormat in field ")(myFieldName)(", failed with:")(realE)("\n");
                res=false;
            }
            
        }
        return true;
    }
}

// object that keeps the journal of the computations done
class PosELogger(T):EmptyObserver!(T){
    PosELoggerGen input;
    LocalSilosI!(T) silos;
    char[] journalPath;
    OutStreamI outStream;
    RLock serialLock;
    string name(){
        return "PosELogger_"~input.myFieldName;
    }
    
    /// creates a new journal logger
    this(PosELoggerGen c,LocalSilosI!(T) silos){
        this.input=c;
        this.silos=silos;
        this.serialLock=new RLock();
        this.outStream=silos.outfileForName(input.logFile,WriteMode.WriteAppend,StreamOptions.CharBase);
    }
    
    /// adds energy for a local point and bCasts addEnergyEval
    void addEnergyEvalLocal(SKey key,Point p,Real energy,Real energyError){
        try{
            this.serialLock.lock();
            scope(exit){ this.serialLock.unlock(); }
            auto s=dumper(&outStream.rawWriteStr);
            s(key)("\t ")(p.data);
            MainPointI!(T) mp=silos.mainPointL(p);
            foreach(x;mp.pos.dynVars.x.sLoop){
                s("\t ")(x,input.posFormat);
            }
            s("\t ")(energy,input.eFormat);
            s("\n");
            if (input.flushEachLine) outStream.flush();
        } catch (Exception e){
            silos.logMsg(delegate void(CharSink s){
                dumper(s)("Warning: Exception while trying to log position and energy from ")(name)(":")(e);
            });
        }
    }
    /// increases the runLevel on all silos, i.e. you should call it only with SKeyVal.All
    /// (at the moment there is no support for dynamic adding/removal of silos)
    void increaseRunLevel(SKey s,RunLevel rLevel){
        outStream.flush();
        if (rLevel>RunLevel.Stopping) outStream.close();
    }
}
