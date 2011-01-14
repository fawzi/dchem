module dchem.pnet.TrackSpecialPoints;
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
import blip.io.StreamConverters;
import tango.io.device.File;
import tango.io.device.Conduit;
import blip.container.Set;
import blip.container.GrowableArray;
import blip.util.NotificationCenter;
import blip.core.Variant;
import dchem.pnet.MainPoint;
import blip.core.stacktrace.StackTrace;

/// a loader of points (what is created by the input)
class TrackSpecialPointsGen:SilosWorkerGen{
    string logfileBaseName="specialPoints";
    bool flushEachLine=false;
    
    mixin(serializeSome("dchem.TrackSpecialPoints",`
    logfileBaseName: base path used for the file where the special points are logged, if emty no log is written (defaults to log)
    flushEachLine: if each line should be flushed`));
    mixin printOut!();
    mixin myFieldMixin!();
    
    this(){}
    
    bool verify(CharSink s){
        return true;
    }
    SilosWorkerI!(Real) silosWorkerReal(){
        auto res=new TrackSpecialPoints!(Real)(this);
        return res;
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        auto res=new TrackSpecialPoints!(LowP)(this);
        return res;
    }
}

class TrackSpecialPoints(T):SilosComponentI!(T){
    Set!(Point) specialPoints;
    LocalSilosI!(T) silos;
    TrackSpecialPointsGen input;
    Callback *myCallback;
    OutStreamI stream;
    
    this(TrackSpecialPointsGen input){
        this.input=input;
        this.specialPoints=new Set!(Point)();
    }
    
    void workOn(LocalSilosI!(T)silos){
        this.silos=silos;
        if (input.logfileBaseName.length>0){
            auto f=new File(input.logfileBaseName~"-"~silos.name~".spLog",File.WriteAppending);
            this.stream=strStreamSyncT!(char)(f);
        }
        myCallback=silos.nCenter.registerCallback("localPointChangedGFlags",
            &flagsChanged,Callback.Flags.ReceiveAll);
        silos.addComponent(this);
    }
    
    char[] name(){
        return "TrackSpecialPoints_"~input.myFieldName;
    }
    char[] kind(){
        return "TrackSpecialPoints";
    }
    
    /// callback that checks the flags changes
    void flagsChanged(cstring notificationName,Callback* callback,Variant oldF){
        assert(notificationName=="localPointChangedGFlags");
        auto flagChange=oldF.get!(GFlagsChange*)();
        auto oldProb=specialPointForGFlags(flagChange.oldGFlags);
        auto newProb=specialPointForGFlags(flagChange.newGFlags);
        if (oldProb != newProb){
            switch(oldProb){
            case Prob.Unlikely,Prob.Possible:
                if (newProb>=Prob.Likely){
                    addSpecialPoint(*flagChange);
                }
                break;
            case Prob.Likely,Prob.Confirmed:
                if (newProb<Prob.Likely){
                    rmSpecialPoint(*flagChange);
                } else {
                    auto oldType=pointTypeForGFlags(flagChange.oldGFlags);
                    auto newType=pointTypeForGFlags(flagChange.newGFlags);
                    if (oldType!=newType){
                        specialPointTypeChange(*flagChange);
                    } else {
                        specialPointProbChange(*flagChange);
                    }
                }
                break;
            default:
                assert(0,"unexpected probability");
            }
        }
    }
    
    void logChange(string change,GFlagsChange flagChange){
        if (stream!is null){
            sinkTogether(&stream.rawWriteStrC,delegate void(CharSink s){
                auto newT=pointTypeForGFlags(flagChange.newGFlags);
                dumper(s)(change)("\t ")(flagChange.point.data)("\t ")(pointTypeStr(newT))("\t ")(newT)("\t ")(specialPointForGFlags(flagChange.newGFlags))("\t ")(flagChange.newGFlags)("\n");
            });
            if (input.flushEachLine) stream.flush();
        }
    }

    void addSpecialPoint(GFlagsChange flagChange){
        bool skip=false;
        synchronized(specialPoints){
            if (specialPoints.contains(flagChange.point)){
                skip=true;
            } else {
                specialPoints.add(flagChange.point);
            }
        }
        if (!skip){
            silos.nCenter.notify("specialPointAdd",Variant(&flagChange));
            logChange("add",flagChange);
        }
    }
    void rmSpecialPoint(GFlagsChange flagChange){
        bool wasRemoved;
        synchronized(specialPoints){
            wasRemoved=specialPoints.remove(flagChange.point);
        }
        if (wasRemoved){
            silos.nCenter.notify("specialPointRm",Variant(&flagChange));
            logChange("rm ",flagChange);
        }
    }
    void specialPointTypeChange(GFlagsChange flagChange){
        silos.nCenter.notify("specialPointTypeChange",Variant(&flagChange));
        logChange("tch",flagChange);
    }
    void specialPointProbChange(GFlagsChange flagChange){
        silos.nCenter.notify("specialPointProbChange",Variant(&flagChange));
        logChange("pch",flagChange);
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
            silos.rmComponentNamed(this.name);
        }
        if (stream!is null){
            stream.close();
            stream=null;
        }
    }
}