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
import blip.container.GrowableArray;
import blip.io.BasicIO;
import dchem.calculator.ProcContext;
import blip.serialization.Serialization;

void clearEF(T)(ParticleSys!(T)pSys){
    pSys.potentialEnergy=T.init;
    pSys.potentialEnergyError=T.init;
    pSys.mddxError=T.init;
    if (pSys.dynVars.mddx.pos !is null)
        pSys.dynVars.mddx.pos[]=T.init;
    // clear also other forces??? for now assume no
}

class PoorManExecuter:CmdTemplateExecuter{
    string energyFile="dchem.energy";
    string forceFile="dchem.forces";

    override CalculationContext getCalculator(bool wait,ubyte[]history){
        auto ctx=new PoorManContext(this,collectAppender(delegate void(CharSink s){
            dumper(s)("PMCtx")(ProcContext.instance.id)("-")(ProcContext.instance.localId.next());
        }));
        return ctx;
    }
    mixin(serializeSome("dchem.PoorManExecuter","executer that can be used to connect any prgram that writes its output in the correct format",`
    energyFile: file where to read energy and energy error (just those numbers in Hartrees, separated by spaces)
    forceFile: file where to read the error on the forces and then the forces in the input order (just numbers in atomic units separated by whitespace)`));
}

/// a context that uses the "poor man" approach: any method works, as long as it puts energy and forces in the appropriate
/// files in the appropriated format (enegy,error),(error,forces)
class PoorManContext:ExecuterContext{
    PoorManExecuter inputPM;
    
    this(PoorManExecuter input, char[] contextId){
        inputPM=input;
        super(input,contextId);
    }
    void readFomattedWhitespaceF(T)(ParticleSys!(T)pSys){
        scope inF=toReaderChar(templateH.targetDir.file(inputPM.forceFile).input);
        scope(exit){ inF.shutdownInput(); }
        scope p=new TextParser!(char)(inF);

        p(pSys.dynVars.mddxError);
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
            throw new Exception("force file '"~inputPM.forceFile~"' is supposed to contain just the error on the forces and the forces as sequence of numbers and nothing else, after reading all forces found '"~tok~"'.",__FILE__,__LINE__);
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
            scope inF=toReaderChar(templateH.targetDir.file(inputPM.energyFile).input);
            scope(exit){ inF.shutdownInput(); }
            scope p=new TextParser!(char)(inF);
            Real ePot,ePotErr;
            p(ePot)(ePotErr);
            if (p.nextToken().length>0){
                throw new Exception("Energy file '"~inputPM.energyFile~" is supposed to contain only energy and its error.",__FILE__,__LINE__);
            }
            switch(activePrecision) {
            case Precision.Real:
                pSysReal.dynVars.potentialEnergy=ePot;
                pSysReal.dynVars.potentialEnergyError=ePotErr;
                break;
            case Precision.LowP:
                pSysLowP.dynVars.potentialEnergy=ePot;
                pSysLowP.dynVars.potentialEnergyError=ePotErr;
                break;
            default:
                assert(0);
            }
        }
    }
}

