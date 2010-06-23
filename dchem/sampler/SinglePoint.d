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

class SinglePoint:Sampler{
    bool calcE;
    bool calcF;
    InputField method;
    mixin myFieldMixin!();

    mixin(serializeSome("dchem.SinglePoint",
    `calcE: the energy should be calculated
    calcF: the forces should be calculated
    method: the method to use`));
    mixin printOut!();
    
    /// does the single point calculation
    void run(){
        auto m=method.method; // method.method;
        if (m is null) throw new Exception("invalid method in field "~myFieldName,__FILE__,__LINE__);
        sout("preparing to calculate\n");
        CalculationContext c=m.getCalculator(true,null);
        sout("will calculate\n");
        c.updateEF(calcE,calcF);
        sout("did calculate\n");
        if (calcE){
            sout("potentialEnergy:")(c.potentialEnergy())("\n");
        }
        if (calcF){
            mixin(withPSys(`
            if (pSys.dynVars.x.dof.length>0 || pSys.dynVars.x.orient.length>0){
                auto w=pSysWriter(pSys);
                sout(w);
            } else {
                auto f=c.mddpos;
                WriteOut.writeXyz!(Real)(sout.call,c.sysStruct,f,"forces "~c.sysStruct.name);
            }`,"c."));
        }
        sout("End calculation\n");
    }
    
    /// stops the calculation if possible
    void stop(){
    }
    
    bool verify(CharSink log){
        bool res=true;
        auto s=dumper(log);
        if (method is null || method.typeId!=InputField.TypeId.Method){
            s("Error: method has to be valid and contain a method")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
}
