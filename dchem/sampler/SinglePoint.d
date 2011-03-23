module dchem.sampler.SinglePoint;
import dchem.Common;
import dchem.sys.ParticleSys;
import dchem.input.RootInput;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import blip.io.Console;
import WriteOut=dchem.input.WriteOut;
import dchem.calculator.FileCalculator;
import dchem.calculator.Calculator;
import dchem.input.WriteOut;
import dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import blip.parallel.mpi.MpiModels;

class SinglePoint:Sampler{
    bool calcE;
    bool calcF;
    InputField method;
    mixin myFieldMixin!();

    mixin(serializeSome("SinglePoint",`Sampler that performs a single point evaluation.`,
    `calcE: the energy should be calculated
    calcF: the forces should be calculated
    method: the method to use`));
    mixin printOut!();
    
    /// does the single point calculation
    void run(LinearComm pWorld,CharSink log){
        auto m=cast(Method)method.contentObj;
        if (m is null) throw new Exception("invalid method in field "~myFieldName,__FILE__,__LINE__);
        sout("method setup");
        m.setup(pWorld,log);
        sout("preparing to calculate\n");
        CalculationContext c=m.getCalculator(true,null);
        sout("will calculate\n");
        c.updateEF(calcE,calcF);
        sout("did calculate\n");
        if (calcE){
            sout("potentialEnergy:")(c.potentialEnergy())("\n");
        }
        auto sysStruct=c.sysStruct;
        if (calcF){
            auto wR0=c.pSysWriterReal();
            auto d=&wR0.mddx.pos.toSegArr!(Real);
            switch(c.activePrecision){
            case Precision.Real:
                auto wR=c.pSysWriterReal();
                if (wR.mddx.dof.data.length>0 || wR.mddx.orient.data.length>0){
                    sout(wR);
                } else {
                    auto fR=wR.mddx.pos.toSegArrFromSysStruct!(Vector!(Real,3))(sysStruct,true);
                    WriteOut.writeXyz!(Real)(sout.call,sysStruct,fR,"forces "~sysStruct.name);
                }
                break;
            case Precision.LowP:
                auto wL=c.pSysWriterLowP();
                if (wL.mddx.dof.data.length>0 || wL.mddx.orient.data.length>0){
                    sout(wL);
                } else {
                    auto fL=wL.mddx.pos.toSegArrFromSysStruct!(Vector!(LowP,3))(sysStruct,true);
                    WriteOut.writeXyz!(LowP)(sout.call,sysStruct,fL,"forces "~c.sysStruct.name);
                }
                break;
            default:
                assert(0);
            }
        }
        sout("End calculation\n");
    }
    
    /// stops the calculation if possible
    void stop(){
    }
    
    bool verify(CharSink log){
        bool res=true;
        auto s=dumper(log);
        if (method is null || (cast(Method)method.contentObj) is null){
            s("Error: method has to be valid and contain a method in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
}
