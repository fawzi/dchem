module dchem.pnet.WaitExplorer;
import dchem.Common;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.io.EventWatcher;
import blip.container.GrowableArray;
import dchem.pnet.EmptyExplorer;

/// defines an explorer that waits the given time before stopping
/// useful just to keep the silos busy
class WaitExplorerDef:SilosWorkerGen{
    double time=10;
    this(){
    }
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.WaitExplorer",`
    time : minutes to wait (10)`));
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        if (time<=0) {
            res=false;
            dumper(s)("time has to be positive in field ")(myFieldName)("\n");
        }
        return res;
    }
    SilosWorkerI!(Real) silosWorkerReal(){
        auto res=new WaitExplorer!(Real)(this);
        return res;
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        auto res=new WaitExplorer!(LowP)(this);
        return res;
    }
}

/// waits the given time before stopping
class WaitExplorer(T):EmptyExplorer!(T){
    void delegate(ExplorerI!(T)) available;
    ev_tstamp endTime;
    GenericWatcher timeWatcher;
    WaitExplorerDef input;
    
    override string name(){
        return "WaitExplorer_"~input.myFieldName;
    }
    
    /// called when the time is up
    void timeIsUp(bool t){
        void delegate(ExplorerI!(T)) mAvailable;
        synchronized(this){
            mAvailable=available;
            available=null;
        }
        if (mAvailable!is null){
            mAvailable(this);
        }
    }
    
    this(WaitExplorerDef input){
        this.input=input;
        this.endTime=cast(ev_tstamp)(ev_time()+input.time*60);
        this.timeWatcher=GenericWatcher.timerCreate(this.endTime,0);
        noToutWatcher.addWatcher(this.timeWatcher,&timeIsUp);
    }
    
    // ExplorerI(T)
    
    /// returns a point to evaluate
    Point pointToEvaluateLocal(SKey k,void delegate(ExplorerI!(T))available,int req){
        if (ev_time()<endTime){
            synchronized(this){
                this.available=available;
            }
            return Point(1);
        }
        return Point(0);
    }
}
