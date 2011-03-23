/// sampler used for debugging purposes
module dchem.sampler.TestSampler;
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
import dchem.util.Rotate;
import blip.narray.NArrayLinAlg;

class TestSampler:Sampler{
    InputField method;
    mixin myFieldMixin!();

    mixin(serializeSome("TestSampler","A sampler that just performs some tests",
    `method: the method to use`));
    mixin printOut!();
    
    /// does the single point calculation
    void run(LinearComm pWorld,CharSink log){
        auto m=cast(Method)method.contentObj;
        if (m is null) throw new Exception("invalid method in field "~myFieldName,__FILE__,__LINE__);
        sout("method setup");
        m.setup(pWorld,log);
        sout("preparing to calculate\n");
        LocalCalculationContext c=cast(LocalCalculationContext)m.getCalculator(true,null);
        // sout("will calculate\n");
        // c.updateEF(calcE,calcF);
        // sout("did calculate\n");
        // if (calcE){
        //     sout("potentialEnergy:")(c.potentialEnergy())("\n");
        // }
        switch(c.activePrecision){
        case Precision.Real:
            doOp(c.pSysReal);
            break;
        case Precision.LowP:
            doOp(c.pSysLowP);
            break;
        default:
                assert(0);
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

    void doOp(T)(ParticleSys!(T) pSys){
        auto d1=pSys.dynVars.dVarStruct.emptyDx();
        d1[]=0;
        d1.pos[LocalPIndex(0,0),0]=Vector!(T,3)(1,0,0);
        auto d2=pSys.dynVars.dVarStruct.emptyDx();
        d2[]=0;
        d2.pos[LocalPIndex(0,0),0]=Vector!(T,3)(0,1,0);
        auto d3=pSys.dynVars.dVarStruct.emptyDx();
        d3[]=0;
        d3.pos[LocalPIndex(0,0),0]=Vector!(T,3)(1,0,0);
        d1=rotateEiV(0,d2,d1);
        sout("after rotation:")(dynPVectorWriter(d1))("\n");
        bypax(d3,d1,cast(T)-1,cast(T)1);
        sout("diff:")(d3.norm2())("\n");
    }
}
