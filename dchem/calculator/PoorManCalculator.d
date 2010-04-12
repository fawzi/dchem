module dchem.calculator.PoorManCalculator;
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

void clearEF(T)(ParticleSys!(T)pSys){
    pSys.potentialEnergy=T.init;
    if (pSys.dynVars.mddx.pos !is null)
        pSys.dynVars.mddx.pos[]=T.init;
    // clear also other forces??? for now assume no
}


/// a context that uses the "poor man" approach: any method works, as long as it puts energy and forces in the appropriate
/// files in the appropriated format
class PoorManContext:ExecuterContext{
    char[] energyFile="dchem.energy";
    char[] forceFile="dchem.forces";
    
    static CalculationContext contextAllocator(CalculationInstance cInstance,Method method,char[]className,char[] contextId){
        sout("pippo1\n");
        auto m=cast(TemplateExecuter)cast(Object)method;
        sout("pippo2\n");
        if (m is null){
            throw new Exception("method must be valid and must be a TemplateExecuter",__FILE__,__LINE__);
        }
        sout("pippo3\n");
        return new PoorManContext(m,cInstance,className,contextId);
    }
    
    this(TemplateExecuter input,CalculationInstance cInstance,char[] className,char[] contextId){
        sout("pippo4\n");
        super(input,cInstance,className,contextId);
        sout("pippo5\n");
    }
    void readFomattedWhitespaceF(T)(ParticleSys!(T)pSys){
        sout("pippo r1\n");
        scope inF=toReaderChar(templateH.targetDir.file(forceFile).input);
        sout("pippo r2\n");
        scope(exit){ sout("will shutdownInput\n"); inF.shutdownInput(); sout("did shutdownInput\n"); }
        sout("pippo r3\n");
        scope p=new TextParser!(char)(inF);

        sout("pippo r4\n");
        auto externalOrder=pSys.sysStruct.externalOrder;
        sout("pippo r5\n");
        pSys.checkMddx();
        sout("pippo r6\n");
        auto f=pSys.dynVars.mddx.pos;
        sout("pippo r7\n");
        foreach (idx;externalOrder.lSortedPIndex.sLoop){
            sout("pippo r8.1\n");
            Vector!(T,3) pos;
            p(pos.x)(pos.y)(pos.z);
            sout("pippo r8.2\n");
            f[LocalPIndex(idx),0]=pos;
        }
        sout("pippo r9\n");
        auto tok=p.nextToken();
        sout("pippo r10\n");
        if (tok.length>0){
            throw new Exception("force file '"~forceFile~"' is supposed to contain just the forces as sequence of numbers and nothing else, after reading all forces found '"~tok~"'.",__FILE__,__LINE__);
        }
        sout("pippo r11\n");
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
    MethodAllocators.defaultAllocators["poorman"]=fun2Dlg(&PoorManContext.contextAllocator);
}
