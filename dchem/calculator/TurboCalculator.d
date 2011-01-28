module dchem.calculator.TurboCalculator;
import dchem.calculator.FileCalculator;
import dchem.calculator.Calculator;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import blip.io.StreamConverters;
import blip.text.TextParser;
import dchem.sys.SegmentedArray;
import dchem.Common;
import dchem.input.RootInput;
import blip.util.TemplateFu:fun2Dlg;
import blip.io.Console;
/+
void clearEF(T)(ParticleSys!(T)pSys){
    pSys.potentialEnergy=T.init;
    pSys.dynVars.mddx[]=0;
}


/// a context that calculates using turbomole
class TurboCalculator:ExecuterContext{
    char[] energyFile="dchem.energy";
    char[] forceFile="dchem.forces";
    
    static CalculationContext contextAllocator(CalculationInstance cInstance,Method method,char[]className,char[] contextId){
        auto m=cast(TemplateExecuter)cast(Object)method;
        if (m is null){
            throw new Exception("method must be valid and must be a TemplateExecuter",__FILE__,__LINE__);
        }
        return new TurboCalculator(m,cInstance,className,contextId);
    }
    
    this(TemplateExecuter input,CalculationInstance cInstance,char[] className,char[] contextId){
        super(input,cInstance,className,contextId);
    }
    void readFomattedWhitespaceF(T)(ParticleSys!(T)pSys){
        scope inF=toReaderChar(templateH.targetDir.file(forceFile).input);
        scope(exit){ inF.shutdownInput(); }
        scope p=new TextParser!(char)(inF);

        auto externalOrder=pSys.sysStruct.externalOrder;
        pSys.checkMddx();
        auto f=pSys.dynVars.mddx.pos;
        foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
            Vector!(T,3) pos;
            p(pos.x)(pos.y)(pos.z);
            f[idx,0]=pos;
        }
        auto tok=p.nextToken();
        if (tok.length>0){
            throw new Exception("force file '"~forceFile~"' is supposed to contain just the forces as sequence of numbers and nothing else, after reading all forces found '"~tok~"'.",__FILE__,__LINE__);
        }
    }
    
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        if (updateF){
            mixin(withPSys("readFomattedWhitespaceF(pSys);"));
        } else {
            if (pSysReal!is null) pSysReal.dynVars.potentialEnergy=Real.init;
            if (pSysLowP !is null) pSysLowP.dynVars.potentialEnergy=LowP.init;
        }
        if (updateE||updateF){
            scope inF=toReaderChar(templateH.targetDir.file(energyFile).input);
            scope(exit){ inF.shutdownInput(); }
            scope p=new TextParser!(char)(inF);
            Real ePot;
            p(ePot);
            if (p.nextToken().length>0){
                throw new Exception("Energy file '"~energyFile~" is supposed to contain only a single number.",__FILE__,__LINE__);
            }
            if (pSysReal!is null) pSysReal.dynVars.potentialEnergy=ePot;
            else if (pSysLowP!is null) pSysLowP.dynVars.potentialEnergy=ePot;
            else assert(0);
        }
    }
}

static this(){
    MethodAllocators.defaultAllocators["turbo"]=fun2Dlg(&TurboCalculator.contextAllocator);
}
+/