module dchem.calculator.TurboCalculator;
import dchem.calculator.Calculator;
import dchem.calculator.FileCalculator;
/+
class TmoleContext:ExecuterContext{
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        assert(0,"to implement in subclasses");
    }
}

static this(){
    MethodAllocators.defaultAllocators["tmole"]=delegate CalculationInstance(CalculationInstance cInstance,Method method,char[] contextId){
        return new TmoleContext(cInstance,method,contextId);
    };
}+/
