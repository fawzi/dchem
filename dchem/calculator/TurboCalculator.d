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
import blip.io.BasicIO;
import blip.container.GrowableArray;
import blip.serialization.Serialization;
import dchem.calculator.ProcContext;
import blip.math.IEEE;
import blip.math.Math;

void clearEF(T)(ParticleSys!(T)pSys){
    pSys.potentialEnergy=T.init;
    pSys.potentialEnergyError=T.init;
    pSys.mddxError=T.init;
    if (pSys.dynVars.mddx.pos !is null)
        pSys.dynVars.mddx.pos[]=T.init;
    // clear also other forces??? for now assume no
}

/// reads a turbomole gradient file in a pSys
struct GradInPSys{
    TextParser!(char)p;
    enum Op{
        Load,  /// values should be loaded in the pSys
        Check, /// values should be checked against the ones in pSys
        Ignore, /// values should be ignored
    }
    Op energyOp; /// operation to perform on the read energy
    Op posOp;    /// operation to perform on the read positons
    Op gradOp;   /// operation to perform on the read gradients
    /// if the pSys is just a contiguos subset starting at the first index up to the second index,
    /// with the third index that gives the total number of particles. Second and third index might
    /// be automatically calculated
    long[] subset;
    /// if an exception should be thrown if element types don't match
    bool throwOnElMismatch=true;
    /// maximum difference in energy (if energyOp==Op.Check)
    Real maxEDiff=1.e-4;
    /// maximum difference in position (if posOp==Op.Check)
    Real maxPosDiff=1.e-3;
    /// maximum difference in gradient (if gradOp==Op.Check)
    Real maxGradDiff=1.e-4;
    /// energy is not compared if the old value is NaN
    bool weakECheck=true;

    /// reads the start of a gradient file (i.e. $grad...)
    void readFileStart(){
        p.newlineIsSpace=true;
        p.skipString("$grad");
        p.newlineIsSpace=false;
        p.skipString("cartesian gradients");
        p.skipLines(1);
    }
    
    /// helper error message
    void fileMsg(CharSink s){
        dumper(s)(" in file ");
        p.parserPos(s);
    }

    /// load a configuration in the given pSys 
    ParticleSys!(T) configInPSys(T)(ParticleSys!(T)pSys,long minCycle=0){
        p.skipWhitespace();
        {
            auto tok=p.nextToken();
            if (tok=="$end") return null;
            if (tok!="cycle") {
                throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("unexpected token '")(tok)("' instead of cycle")(&fileMsg);
                        }),__FILE__,__LINE__);
            }
        }
        p.skipString("=");
        long cycle;
        p(cycle);
        while (cycle<minCycle){
            if (!p.skipLines(1,false)){
                throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("found EOF while trying to find cycle ")(minCycle)(&fileMsg);
                        }),__FILE__,__LINE__);
                
            }
            p.skipWhitespace();
            auto token=p.nextToken();
            if (token=="$end"){
                throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("found $end while trying to find cycle ")(minCycle)(&fileMsg);
                        }),__FILE__,__LINE__);
                
            } else if (token=="cycle"){
                p.skipString("=");
                p(cycle);
            }
        }
        pSys.iteration=cycle;
        if (!((energyOp==Op.Ignore) || (weakECheck && energyOp==Op.Check && isNaN(pSys.dynVars.potentialEnergy)))){
            p.skipString("SCF energy =");
            Real eNow;
            p(eNow);
            switch(energyOp){
            case Op.Check:
                if (!(abs(pSys.dynVars.potentialEnergy-eNow)<maxEDiff)){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("energy change was too big:")(abs(pSys.dynVars.potentialEnergy-eNow))
                                    (" (")(eNow)(" vs ")(pSys.dynVars.potentialEnergy)(")")(&fileMsg);
                            }),__FILE__,__LINE__);
                }
                break;
            case Op.Load:
                pSys.dynVars.potentialEnergy=eNow;
                break;
            case Op.Ignore:
                break;
            default:
                assert(0);
            }
        }
        long iParticle=0;
        if (subset.length>0){
            for (long i=0;i<subset[0];++i){
                p.skipLines(1);
                ++iParticle;
                T x;
                p(x);
            }
        }
        auto externalOrder=pSys.sysStruct.externalOrder;
        if (subset.length>1 && subset[1]-subset[0]
            !=externalOrder.nLocalParticles(pSys.sysStruct.levels[0])){
            throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("number of particles in subset is different from the ones in particleSystem")
                            (subset[1]-subset[0])(" vs ")(externalOrder.nLocalParticles(pSys.sysStruct.levels[0]));
                            fileMsg(s);
                    }),__FILE__,__LINE__);
        }
        switch(posOp){
        case Op.Load:
            pSys.checkX();
            // pass to the next...
        case Op.Check:
            SegmentedArray!(Vector!(T,3)) xPos=pSys.dynVars.x.pos;
            auto maxPosDiff2=pow2(maxPosDiff);
            foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
                p.skipLines(1);
                ++iParticle;
                char[] symb;
                Vector!(T,3) pos;
                p(pos.x)(pos.y)(pos.z);
                if (posOp==Op.Check){
                    if (xPos!is null && (!((pos-xPos[idx,0]).norm22<maxPosDiff2))){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                                    dumper(s)("difference in position of particle ")(iParticle)
                                        (" is too large ")(pos)(" vs ")(xPos[idx,0])(&fileMsg);
                                }),__FILE__,__LINE__);
                    }
                } else {
                    xPos[idx,0]=pos;
                }
                p(symb,false);
                if (throwOnElMismatch){
                    if (symb.length>0 && symb[0]>='a'&&symb[0]<='z')
                        symb[0]=symb[0]+'A'-'a';
                    if (symb.length>1 && symb[1]>='A'&&symb[1]<='Z')
                        symb[1]=symb[1]-'A'+'a';
                    if (symb!=pSys.sysStruct.kinds[cast(size_t)idx.kind].symbol){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                                    dumper(s)("inconsistent Symbol detected:")(symb)(" vs ")
                                        (pSys.sysStruct.kinds[cast(size_t)idx.kind].symbol)(&fileMsg);
                                }),__FILE__,__LINE__);
                    }
                }
            }
            break;
        case Op.Ignore:
            foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
                p.skipLines(1);
                ++iParticle;
                char[] symb;
                Vector!(T,3) pos;
                p(pos.x)(pos.y)(pos.z)(symb);
            }
            break;
        default:
            assert(0);
        }
        // reads until grad starts, handle also the case where there are no particles???
        Vector!(T,3) pos;
        try{
            char[]tok;
            do{
                ++iParticle;
                p.skipLines(1);
                p(pos.x)(pos.y)(pos.z);
                p.skipWhitespace();
                tok=p.nextToken();
            } while (tok.length!=0);
            --iParticle;
        } catch(Exception e){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("error trying to find the start of the gradient part ")(p);
            }),__FILE__,__LINE__,e);
        }
        long iParticle2=0;
        if (subset.length!=0){
            for (long i=0;i<subset[0];++i){
                ++iParticle2;
                p.skipLines(1);
                p(pos.x)(pos.y)(pos.z);
            }
        }
        // read the gradients...
        if (gradOp==Op.Load) pSys.checkMddx();
        auto maxFDiff2=pow2(maxGradDiff);
        SegmentedArray!(Vector!(T,3)) f=pSys.dynVars.mddx.pos;
        bool readPos=false;
        foreach (idx;externalOrder.gSortedLocalPIndex.sLoop){
            if (readPos) p(pos.x)(pos.y)(pos.z);
            readPos=true;
            ++iParticle2;
            if (gradOp==Op.Load) {
                f[idx,0]=-pos;
            } else if (gradOp==Op.Check){
                if (f!is null && (!((pos+f[idx,0]).norm22<maxFDiff2))){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("difference in force of particle ")(iParticle)
                                    (" is too large ")((pos+f[idx,0]).norm22)(" (")
                                    (-pos)(" vs ")(f[idx,0])(")")(&fileMsg);
                            }),__FILE__,__LINE__);
                }
            }
            p.skipLines(1);
        }
        while (iParticle2<iParticle){
            ++iParticle2;
            if (readPos) p(pos.x);
            readPos=true;
            p.skipLines(1);
        }
        p.skipWhitespace();
        auto tok=p.peek(&p.scanString);
        if (tok!="$end" && tok!="cycle"){
            throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("after reading a configuration found '")(tok)("'")(&fileMsg);
                    }),__FILE__,__LINE__);
        }
        return pSys;
    }
    
    static GradInPSys opCall(TextParser!(char)p){
        GradInPSys res;
        res.p=p;
        p.newlineIsSpace=false;
        return res;
    }

    /// deallocates the data that might have been allocated by this
    void deallocData(){
    }
}

/// executer that uses turbomole to calculate the energy and forces
class TurboExecuter:CmdTemplateExecuter{
    string enerFile="energy";
    string eErrorFile="ridft.out";
    string gradFile="gradient";
    long lastCycle=0;
    long[] subset;
    Real maxEDiff=1.e-4;
    Real maxPosDiff=1.e-3;
    bool incrementCycle=true;
    bool throwOnElMismatch=true;
    bool shouldBeLast=true;
    bool executeF0DoesNotIncrement=true;

    override CalculationContext getCalculator(bool wait,ubyte[]history){
        auto ctx=new TurboContext(this,collectAppender(delegate void(CharSink s){
            dumper(s)("TMCtx")(ProcContext.instance.id)("-")(ProcContext.instance.localId.next());
        }));
        return ctx;
    }
    mixin(serializeSome("dchem.TurboExecuter",`
        enerFile: file where the energy is stored (energy)
        eErrorFile: file where to read the convergence error (ridft.out)
        gradFile: file where to read the gradient(gradient)
        lastCycle: last cycle for energy and gradient (0)
        incrementCycle: if the cycle should be incremented
        subset: if just a subset [from..to) 0 based of the atoms should be considered
        maxEDiff: maximum difference in energy in au (1.e-4)
        maxPosDiff: maximum difference in the positions in au (1.e-3)
        throwOnElMismatch: if an exception should be raised on a mismatch in the element (atom) symbol (true)
        shouldBeLast: if the gradeint/energy read should always be the last (true)
        executeF0DoesNotIncrement: if executeF0 (if given) does not increment the cycle (true)`));
    override bool verify(CharSink s){
        bool res=super.verify(s);
        if (subset.length!=0 && subset.length!=2){
            dumper(s)("subset in field ")(myFieldName)(" must be either empty, or two number representing a 0 based range [startInclusive,endExclusive)\n");
            res=false;
        }
        return res;
    }

}

/// a context that uses the energy and gradients in the turbomole format
class TurboContext:ExecuterContext{
    TurboExecuter inputTurbo;
    long expectedEPos;
    long expectedGPos;

    this(TurboExecuter input, char[] contextId){
        inputTurbo=input;
        super(input,contextId);
    }
    
    void readGradient(T)(ParticleSys!(T)pSys){
        try{
            scope inF=toReaderChar(templateH.targetDir.file(inputTurbo.gradFile).input);
            scope(exit){ inF.shutdownInput(); }
            scope p=new TextParser!(char)(inF);
            auto gParser=GradInPSys(p);
            gParser.energyOp=GradInPSys.Op.Check;
            gParser.posOp=GradInPSys.Op.Check;
            gParser.gradOp=GradInPSys.Op.Load;
            gParser.subset=inputTurbo.subset;
            gParser.maxEDiff=inputTurbo.maxEDiff;
            gParser.maxPosDiff=inputTurbo.maxPosDiff;
            gParser.throwOnElMismatch=inputTurbo.throwOnElMismatch;
            gParser.readFileStart();
            auto pSysRes=gParser.configInPSys(pSys,expectedGPos);
            if (expectedGPos==-1){
                auto pSysNew=pSysRes;
                while (pSysNew!is null){
                    pSysRes=pSysNew;
                    pSysNew=gParser.configInPSys(pSys,expectedGPos);
                }
            } else {
                if (pSysRes is null || pSysRes.iteration!=expectedGPos && expectedGPos!=0){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("could not find cycle ")(expectedGPos)(&gParser.fileMsg);
                            }),__FILE__,__LINE__);
                }
            }
            if (inputTurbo.shouldBeLast && gParser.configInPSys(pSys)!is null){
                throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("cycle ")(expectedGPos)(" is not the last configuration")(&gParser.fileMsg);
                        }),__FILE__,__LINE__);
            }
            gParser.deallocData();
        } catch (Exception e) {
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("Error in ")(&desc)(" trying to read gradient");
            }),__FILE__,__LINE__,e);
        }
    }

    void readEnergy(T)(ParticleSys!(T)pSys){
        try{
            scope inF=toReaderChar(templateH.targetDir.file(inputTurbo.enerFile).input);
            scope(exit){ inF.shutdownInput(); }
            scope p=new TextParser!(char)(inF);
            p.skipString("$energy");
            p.newlineIsSpace=false;
            p.skipLines(1);
            long cycle;
            p(cycle);
            while(cycle<expectedEPos){
                p.skipLines(1);
                p(cycle);
            }
            Real eAtt;
            p(eAtt);
            if (expectedEPos==-1){
                while (true){
                    p.skipLines(1);
                    if (p.skipString("$end",false)) break;
                    p(cycle)(eAtt);
                }
            } else {
                if (cycle!=expectedEPos && expectedEPos!=0){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("could not find cycle ")(expectedEPos)(&p.parserPos);
                            }),__FILE__,__LINE__);
                }
                if (inputTurbo.shouldBeLast){
                    p.skipLines(1);
                    if (!p.skipString("$end",false)){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                                    dumper(s)("cycle ")(expectedEPos)(" is not the last configuration")(&p.parserPos);
                                }),__FILE__,__LINE__);
                    }
                }
            }
            pSys.dynVars.potentialEnergy=eAtt;
            pSys.iteration=cycle;
        } catch (Exception e){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("Error in ")(&desc)(" trying to read energy");
            }),__FILE__,__LINE__,e);
        }
    }

    Real readEError(){
        try{
            if (inputTurbo.eErrorFile.length>0){
                scope inF=toReaderChar(templateH.targetDir.file(inputTurbo.eErrorFile).input);
                scope(exit){ inF.shutdownInput(); }
                scope p=new TextParser!(char)(inF);
                Real errFock,errFiaBlock;
                while(true){
                    if (p.skipString("max. resid. ",false)){
                        if (p.skipString("norm for Fia-block",false)&&p.skipString("=",false))
                            p(errFiaBlock);
                        if (p.skipString("fock norm",false) && p.skipString("=",false))
                            p(errFock);
                    }
                    if (p.skipLines(1,false)==0) break;
                }
                switch((isNaN(errFock)?2:0)+(isNaN(errFiaBlock)?1:0)){
                case 3:
                    return sqrt(pow2(errFock)+pow2(errFiaBlock));
                case 2:
                    return errFock;
                case 1,0:
                    return errFiaBlock;
                default: assert(0);
                }
            }
        } catch (Exception e){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("Error in ")(&desc)(" trying to read gradient");
            }),__FILE__,__LINE__,e);
        }
        return Real.init;
    }
    
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        if (inputTurbo.incrementCycle){
            if(!((!updateE)&&updateF&&changeLevel==ChangeLevel.NoChange&&
                inputTurbo.executeF0.length!=0&&inputTurbo.executeF0DoesNotIncrement))
            {
                ++expectedEPos;
                expectedGPos=expectedEPos;
            }
        }
        if (updateE){
            mixin(withPSys(`
            readEnergy(pSys);
            pSys.dynVars.potentialEnergyError=readEError();
            `));
        }
        if (updateF){
            mixin(withPSys(`
            readGradient(pSys);
            auto eErr=pSys.dynVars.potentialEnergyError;
            if (!isNaN(eErr)){
                pSys.dynVars.mddxError=sqrt(eErr);
            } else {
                pSys.dynVars.mddxError=Real.init;
            }
            `));
        }
    }
}
