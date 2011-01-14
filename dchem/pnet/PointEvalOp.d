/// an eval op that evaluates a single point
module dchem.pnet.PointEvalOp;
import dchem.Common;
import dchem.calculator.CalculatorModels;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.serialization.Serialization;
import blip.io.EventWatcher: ev_tstamp,ev_time;

/// an operation on a context and SilosConnection
class PointEvalOp(T):EvalOp!(T){
    Point point;
    bool evalGrad=false;

    this(){ super(); }
    this(Point point,bool evalGrad=false){
        super();
        this.point=point;
        this.evalGrad=evalGrad;
    }
    /// performs the operation
    void doOp(){
        if (status<Status.Success){
            assert(silos!is null);
            bool success=false;
            try{
                auto time=ev_time();
                auto localP=silos.createLocalPoint(point,time);
                scope(exit) { localP.release(); }
                if (localP.evalWithContext(ctx,evalGrad)){
                    version(TrackWorkAsker) {
                        sinkTogether(sout,delegate void(CharSink s){
                            dumper(s)(this)(" did eval point ")(point)("\n");
                        });
                    }
                } else {
                    version(TrackWorkAsker) {
                        silos.logMsg(delegate void(CharSink s){
                            dumper(s)(this)(" point ")(point)(" was already evaluated\n");
                        });
                    }
                }
                success=true;
            } catch (Exception e){
                silos.logMsg(delegate void(CharSink s){
                    dumper(s)(this)(" evaluation failed with exception ")(e)("\n");
                });
            }
            if (success){
                updateStatus(Status.Success);
                silos.updateEvalStatus(owner,id,null,Status.Success);
            } else {
                updateStatus(Status.Failure);
                silos.updateEvalStatus(owner,id,null,Status.Failure);
                silos.evaluationFailed(SKeyVal.All,point);
            }
        } else {
            silos.logMsg(delegate void(CharSink s){
                dumper(s)("doOp called on EvalOp ")(id)(" which has already run status ")(status);
            });
        }
    }
    /// tries to stop the operation early
    void stopOp(){ }
    /// update the status
    bool updateStatus(Status status){
        synchronized(this){
            if (status==Status.Failure){
                ++attempts;
            }
            if (status>=this.status){
                this.status=status;
            }
        }
        return (attempts<=maxAttempts);
    }
    
    mixin(serializeSome("PointEvalOp!("~T.stringof~")",`point`));
}
