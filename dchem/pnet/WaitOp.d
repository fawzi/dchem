/// an eval op that evaluates a single point
module dchem.pnet.WaitOp;
import dchem.Common;
import dchem.calculator.CalculatorModels;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.serialization.Serialization;
import blip.io.EventWatcher: ev_tstamp,ev_time,noToutWatcher;

/// an operation on a context and SilosConnection
class WaitOp(T):EvalOp!(T){
    double toWait=10;

    this(){ super(); }
    this(double toWait){
        super();
        this.toWait=toWait;
    }
    /// performs the operation
    void doOp(){
        if (status<Status.Success){
            noToutWatcher.sleepTask(toWait);
        }
        updateStatus(Status.Success);
        silos.updateEvalStatus(owner,id,null,Status.Success);
    }
    /// tries to stop the operation early
    void stopOp(){ }
    /// update the status
    bool updateStatus(Status status){
        super.updateStatus(status);
        return false;
    }
    
    mixin(serializeSome("WaitOp!("~T.stringof~")",`An operation that just waits.`,`toWait`));
}
