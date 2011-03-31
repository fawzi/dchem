/// logs position and energy of the points
module dchem.pnet.PrintAttractors;
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
class PrintAttractorsGen:SilosWorkerGen{
    char[] fileName="AttractorLog.log";
    char[] posFormat="f8.4";
    char[] eFormat="g12.6";
    bool printPos=false;
    bool printEnergy=true;
    this(){}
    SilosWorkerI!(Real) silosWorkerReal(){
        return new PrintAttractors!(Real)(this);
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        return new PrintAttractors!(LowP)(this);
    }
    mixin(serializeSome("PrintAttractors",`Prints the attractor of each point`,
    `fileName: base filename used for the logs
    posFormat: format string to format the positions
    eFormat: format string to format the energy
    printPos: if the position of the point should be printed (false)
    printEnergy: if the energy should be printed (true)`));
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

// object that prints the attractors
class PrintAttractors(T):SilosWorkerI!(T){
    PrintAttractorsGen input;
    LocalSilosI!(T) silos;
    char[] journalPath;
    OutStreamI outStream;
    RLock serialLock;
    string name(){
        return "PrintAttractors"~input.myFieldName;
    }
    
    /// creates a new journal logger
    this(PrintAttractorsGen c){
        this.input=c;
        this.serialLock=new RLock();
    }
    
    /// adds energy for a local point and bCasts addEnergyEval
    void workOn(LocalSilosI!(T) silos){
        this.silos=silos;
        this.serialLock.lock();
        this.outStream=silos.outfileForName(input.fileName,WriteMode.WriteAppend,StreamOptions.CharBase);
        scope(exit){ this.serialLock.unlock(); }
        auto s=dumper(&outStream.rawWriteStr);
        foreach (p,mp;silos){
            Attractor attractor;
            synchronized(mp){
                attractor=mp.attractor;
            }
            s(p.data)("\t ")(attractor.throughPoint.data)("\t ")(attractor.minimum.data);
            if (input.printEnergy) s("\t ")(mp.energy,input.eFormat);
            if (input.printPos){
                foreach(x;mp.pos.dynVars.x.sLoop){
                    s("\t ")(x,input.posFormat);
                }
            }
            s("\n");
        }
        outStream.flush();
    }
}
