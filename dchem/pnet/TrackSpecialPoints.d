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
import blip.io.FileStream;
import blip.container.Set;
import blip.container.GrowableArray;
import blip.util.NotificationCenter;
import blip.core.Boxer;
import dchem.pnet.MainPoint;
import blip.core.stacktrace.StackTrace;
import dchem.calculator.FileCalculator;
import dchem.calculator.Calculator;
import blip.core.Array;

/// a loader of points (what is created by the input)
class TrackSpecialPointsGen:SilosWorkerGen{
    string logFileName="specialPoints.log";
    bool flushEachLine=true;
    bool logAllGFlagsChanges=false;
    EvalLog[] configLogs=[{targetFile:"specialPoints.xyz",format:"xyz"}];
    
    mixin(serializeSome("TrackSpecialPoints",`Writes out the special points found so far`,
    `logFileName: base path used for the file where the special points are logged, if emty no log is written (defaults to specialPoints.log)
    flushEachLine: if each line should be flushed (true)
    logAllGFlagsChanges: if all gFlags changes should be logged (default is false)`));
    mixin printOut!();
    mixin myFieldMixin!();
    
    this(){}
    
    bool verify(CharSink s){
        bool res=true;
        foreach(l;configLogs){
            if (l.targetFile.length==0){
                dumper(s)("missing targetFile in configLogs in field ")(myFieldName)("\n");
                res=false;
            }
            if (find(WriteOut.writeConfigFormats,l.format)==WriteOut.writeConfigFormats.length){
                auto w=dumper(s);
                w("format in onELog has to be one of ");
                foreach(i,f;WriteOut.writeConfigFormats){
                    if (i!=0) w(", ");
                    w(f);
                }
                w(" and not '")(l.format)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        return res;
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
        if (input.logFileName.length>0){
            this.stream=silos.outfileForName(input.logFileName,WriteMode.WriteAppend,StreamOptions.CharBase|StreamOptions.Sync);
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
    void flagsChanged(cstring notificationName,Callback* callback,Box oldF){
        assert(notificationName=="localPointChangedGFlags");
        auto flagChange=unbox!(GFlagsChange*)(oldF);
        auto oldProb=specialPointForGFlags(flagChange.oldGFlags);
        auto newProb=specialPointForGFlags(flagChange.newGFlags);
        bool logged=false;
        if (oldProb != newProb){
            switch(oldProb){
            case Prob.Unlikely,Prob.Possible:
                if (newProb>=Prob.Likely){
                    logged=true;
                    addSpecialPoint(*flagChange);
                }
                break;
            case Prob.Likely,Prob.Confirmed:
                logged=true;
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
        } else if (oldProb>=Prob.Likely){
            auto oldType=pointTypeForGFlags(flagChange.oldGFlags);
            auto newType=pointTypeForGFlags(flagChange.newGFlags);
            if (oldType!=newType){
                logged=true;
                specialPointTypeChange(*flagChange);
            }
        }
        if ((!logged) && input.logAllGFlagsChanges){
            logChange("gch",*flagChange);
        }
    }
    
    void logChange(string change,GFlagsChange flagChange){
        if (stream!is null){
            sinkTogether(&stream.rawWriteStrC,delegate void(CharSink s){
                auto newT=pointTypeForGFlags(flagChange.newGFlags);
                dumper(s)(change)("\t ")(flagChange.point.data)("\t ")(pointTypeStr(newT))("\t ")(newT)("\t ")(specialPointForGFlags(flagChange.newGFlags))("\t ")(flagChange.oldGFlags)("\t ")(flagChange.newGFlags)("\n");
            });
            if (input.flushEachLine) stream.flush();
        }
        foreach(l;input.configLogs){
            auto f=silos.outfileForName(l.targetFile,WriteMode.WriteAppend);
            scope(exit){ f.flush(); f.close(); }
            auto lp=silos.mainPointL(flagChange.point);
            auto externalRef=collectAppender(delegate void(CharSink s){
                                                 dumper(s)(flagChange.point.data);
                                             });
            WriteOut.writeConfig(f,lp.pos,l.format,externalRef);
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
            silos.nCenter.notify("specialPointAdd",box(&flagChange));
            logChange("add",flagChange);
        }
    }
    void rmSpecialPoint(GFlagsChange flagChange){
        bool wasRemoved;
        synchronized(specialPoints){
            wasRemoved=specialPoints.remove(flagChange.point);
        }
        if (wasRemoved){
            silos.nCenter.notify("specialPointRm",box(&flagChange));
            logChange("rm ",flagChange);
        }
    }
    void specialPointTypeChange(GFlagsChange flagChange){
        silos.nCenter.notify("specialPointTypeChange",box(&flagChange));
        logChange("tch",flagChange);
    }
    void specialPointProbChange(GFlagsChange flagChange){
        silos.nCenter.notify("specialPointProbChange",box(&flagChange));
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
