module dchem.calculator.Cp2kCalculator;
import dchem.calculator.Calculator;
import dchem.calculator.FileCalculator;
import dchem.sys.ParticleSys;
import dchem.input.RootInput;
/+
class Cp2kContext:ExecuterContext{
    int port;
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        assert(0,"to implement in subclasses");
    }
    this(ClassInstanceManager manager,char[]instanceId,char[]folder,uint maxContexts){
        super(manager,instanceId,folder,maxContexts);
    }
}

static this(){
    MethodAllocators.defaultAllocators["cp2k"]=delegate CalculationInstance(CalculationInstance cInstance,Method method,char[] contextId){
        return new Cp2kContext(cInstance,method,contextId);
    };
}
+/