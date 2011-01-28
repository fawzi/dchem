module dchem.pnet.AddRefPos;
import dchem.Common;
import dchem.input.RootInput;
import dchem.pnet.PNetModels;
import dchem.sys.ParticleSys;
import dchem.pnet.PointEvalOp;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import dchem.pnet.DirArray;

/// adds the ref pos to the silos as point to evaluate
class AddRefPosGen:SilosWorkerGen {
    this(){}
    SilosWorkerI!(Real) silosWorkerReal(){
        return new AddRefPos!(Real)(this);
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        return new AddRefPos!(LowP)(this);
    }
    mixin(serializeSome("dchem.AddRefPos",""));
    mixin printOut!();
    mixin myFieldMixin!();
    bool verify(CharSink s){
        return true;
    }
}

class AddRefPos(T):SilosWorkerI!(T) {
    AddRefPosGen input;
    this(AddRefPosGen input){
        this.input=input;
    }
    
    void workOn(LocalSilosI!(T) silos){
        silos.logMsg1("AddRefPos");
        assert(silos!is null);
        if (silos.paraEnv.myRank==0){
            ParticleSys!(T) pos=silos.refPos();
            assert(pos!is null);
            auto newP=silos.newPointAt(pos.dynVars.x,Point(0));
            assert(newP!is null);
            newP=silos.bcastPoint(newP);
            newP.exploredDirs.atomicCAS(0,DirFlags.Explored,DirFlags.Free);
            EvalOp!(T) newOp=new PointEvalOp!(T)(newP.point,newP.shouldCalculateGradient());
            silos.registerPendingOp(newOp);
            silos.addEvalOp(SKeyVal.Master,newOp,true);
        }
    }
}
