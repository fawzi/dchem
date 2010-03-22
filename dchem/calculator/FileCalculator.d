module dchem.calculator.FileCalculator;
import blip.io.BasicIO;
import blip.io.StreamConverters;
import tango.io.vfs.FileFolder;
import blip.t.time.Time;
import blip.t.time.Clock;
import WriteOut=dchem.input.WriteOut;
import dchem.calculator.Calculator;
import dchem.input.RootInput;
import blip.sync.UniqueNumber;
import blip.serialization.Serialization;
import tango.sys.Process;
import blip.container.GrowableArray;
import Path=tango.io.Path;
import blip.io.Console;

/// represent a manager for local processes
class LocalExeInstanceManager:ClassInstanceManager {
    UniqueNumber!(size_t) lastInstanceId;
    char[] baseDir;
    VfsFolder _baseDirectory;
    bool registered=false;
    
    mixin(serializeSome("dchem.LocalExecuter","baseDir|maxInstances|addToCache"));
    this(){}
    this(char[]baseDir,size_t maxCalc=1,bool addToCache=true){
        super(maxCalc,addToCache);
        this.baseDir=baseDir;
        this.lastInstanceId=UniqueNumber!(size_t)(1);
    }
    VfsFolder baseDirectory(){
        if (_baseDirectory is null){
            _baseDirectory=new FileFolder(baseDir,true);
        }
        return _baseDirectory;
    }
    override CalculationInstance newInstance(uint maxContexts){
        char[] calcInstanceId=collectAppender(delegate void(CharSink s){
            s(myFieldName);
            s("-");
            writeOut(s,lastInstanceId.next());
        });
        auto newF=baseDirectory.folder(calcInstanceId).create();
        char[] dir=baseDir;
        sout("LocalExeInstanceManager will create new Instance pippo\n");
        return new ExecuterInstance(this,calcInstanceId,newF.toString(),maxContexts);
    }
    bool verify(CharSink sink){
        auto s=dumperP(sink);
        if (baseDir.length==0) {
            s("Warning: no baseDir given in field ")(myFieldName)(", assuming '.'\n");
            baseDir=".";
        } else if (! Path.exists(baseDir)){
            s("Warning: non existing baseDir '")(baseDir)("' given in field ")(myFieldName)("\n");
        }
        return true;
    }
}
/// represent one local execution process
class ExecuterInstance:CalcInstance{
    char[] baseDir;
    VfsFolder _baseDirectory;
    Time creationTime;
    
    this(ClassInstanceManager manager,char[]instanceId,char[]folder,uint maxContexts){
        super(manager,instanceId,maxContexts);
        baseDir=folder;
        creationTime=Clock.now;
    }
    VfsFolder baseDirectory(){
        if (_baseDirectory is null){
            _baseDirectory=new FileFolder(baseDir,true);
        }
        return _baseDirectory;
    }
}

class TemplateExecuter: Method {
    InputField superTemplate;
    InputField executer;
    InputField startConfig;
    char[] templateDir;
    char[] setupCmd, executeInitialE, executeStructChangeE, executePosChangeE, executeSmallPosChangeE, executeSmoothPosChangeE, executeDefaultE;
    char[] stopCmd, executeInitialEF, executeStructChangeEF, executePosChangeEF, executeSmallPosChangeEF, executeSmoothPosChangeEF, executeDefaultEF;
    char[][char[]] subs;
    char[] methodClass;
    bool writeReplacementsDict;
    bool overwriteUnchangedPaths;
    int maxContexts=1;
    VfsFolder _templateDirectory;
    
    VfsFolder templateDirectory(){
        if (_templateDirectory is null){
            if (templateDir) _templateDirectory=new FileFolder(templateDir);
        }
        return _templateDirectory;
    }
    bool verify(CharSink sink){
        bool res=true;
        auto s=dumperP(sink);
        if (superTemplate !is null){
            if ((cast(TemplateExecuter)superTemplate.content)is null){
                s("Error: superTemplate should be of type dchem.TemplateExecuter in field "~myFieldName);
                res=false;
            }
        }
        return res;
    }
    
    CalculationInstance getCalcInstance(InstanceGetFlags flags=InstanceGetFlags.ReuseWait){
        sout("entering TemplateExecuter.getCalcInstance\n");
        if (executer is null){
            throw new Exception("Expected a valid executer in the executer field of the dchem.TemplateExecuter "~myFieldName,__FILE__,__LINE__);
        }
        auto manager=cast(ClassInstanceManager)cast(Object)executer.content;
        if (manager is null){
            throw new Exception("Expected an executer (subclass or ClassInstanceManager) in field "~myFieldName,__FILE__,__LINE__);
        }
        sout("call manager.getCalcInstance\n");
        return manager.getInstance(flags,maxContexts);
    }
    CalculationContext getCalculator(CalculationInstance i){
        return i.newContext(this,methodClass,MethodAllocators.defaultAllocators[methodClass]);
    }
    char[] setupCommand(){
        if (setupCmd.length==0 && superTemplate!is null){
            auto te=cast(TemplateExecuter)superTemplate;
            if (te!is null){
                return te.setupCommand();
            }
        }
        return setupCmd;
    }
    char[] stopCommand(){
        if (stopCmd.length==0 && superTemplate!is null){
            auto te=cast(TemplateExecuter)superTemplate;
            if (te!is null){
                return te.stopCommand();
            }
        }
        return stopCmd;
    }
    char[] commandFor(bool energy,bool force,int changeLevel){
        char[] res;
        if (force==false){
            if (energy==false) return "NONE";
            switch (changeLevel){
            case ChangeLevel.FirstTime:
                res=executeInitialE;
                break;
            case ChangeLevel.AllChanged:
                res=executeStructChangeE;
                break;
            case ChangeLevel.PosChanged:
                res=executePosChangeE;
                break;
            case ChangeLevel.SmallPosChange:
                res=executeSmallPosChangeE;
            break;
            case ChangeLevel.SmoothPosChange:
                res=executeSmoothPosChangeE;
            break;
            default:
                throw new Exception(collectAppender(delegate void(CharSink s){ s("unexpected changeLevel "); writeOut(s,changeLevel); }),
                    __FILE__,__LINE__);
            }
            if (res.length==0){
                res=executeDefaultE;
            }
        } else {
            switch (changeLevel){
            case ChangeLevel.FirstTime:
                res=executeInitialEF;
                break;
            case ChangeLevel.AllChanged:
                res=executeStructChangeEF;
                break;
            case ChangeLevel.PosChanged:
                res=executePosChangeEF;
                break;
            case ChangeLevel.SmallPosChange:
                res=executeSmallPosChangeEF;
            break;
            case ChangeLevel.SmoothPosChange:
                res=executeSmoothPosChangeEF;
            break;
            default:
                throw new Exception(collectAppender(delegate void(CharSink s){ s("unexpected changeLevel "); writeOut(s,changeLevel); }),
                    __FILE__,__LINE__);
            }
            if (res.length==0){
                res=executeDefaultE;
            }
        }
        if (res.length==0 && superTemplate !is null){
            auto te=cast(TemplateExecuter)superTemplate.content;
            if (te!is null){
                res=te.commandFor(energy,force,changeLevel);
            }
        }
        return res;
    }
    void addFullSubs(char[][char[]] sub){
        assert(sub!is subs,"cannot add directly to own subs");
        if (superTemplate!is null){
            auto te=cast(TemplateExecuter)superTemplate;
            if (te!is null){
                te.addFullSubs(sub);
            }
        }
        foreach (k,v;subs){
            sub[k]=v;
        }
    }
    // serialization stuff
    mixin(serializeSome("dchem.TemplateExecuter",
    `superTemplate: a template where to take default values
    executer: the executer allocating resources for this method
    methodClass: class of the method to use (i.e. program, needed to know how to get energy,forces,...)
    startConfig: the initial configuration
    templateDir: where to find the definition of the template (tipically a directory)
    subs: keyword and their substitutions to apply to the templates (as dictionary string -> string)
    setupCmd: a command that should be executed to set up the context
    stopCmd: a command that should be executed to stop the context
    executeInitialE: command executed to calculate the first energy, set to NONE to deactivate
    executeStructChangeE: command executed to calculate the energy when the structure is changed, set to NONE to deactivate
    executePosChangeE: command executed to calculate the energy when only positions changed, set to NONE to deactivate
    executeSmallPosChangeE: command executed to calculate the energy when positions changed by a small amount, set to NONE to deactivate
    executeSmoothPosChangeE: command executed to calculate the energy when positions changed smoothly, set to NONE to deactivate
    executeDefaultE: default command for all execute*E commands (if left empty)
    executeInitialEF: command executed to calculate the first energy and forces, set to NONE to deactivate
    executeStructChangeEF: command executed to calculate the energy and forces when the structure is changed, set to NONE to deactivate
    executePosChangeEF: command executed to calculate the energy and forces when only positions changed, set to NONE to deactivate
    executeSmallPosChangeEF: command executed to calculate the energy and forces when positions changed by a small amount, set to NONE to deactivate
    executeSmoothPosChangeEF: command executed to calculate the energy and forces when positions changed smoothly, set to NONE to deactivate
    executeDefaultEF: default command for all execute*EF commands (if they are empty)
    writeReplacementsDict: if a dictionary with the replacements performed should be written (default is false)
    overwriteUnchangedPaths: if paths that are already there should be overwitten (default is false)
    maxContexts: maximum number of contexts per calculation instance (default is 1)
    `));
    mixin printOut!();
    mixin myFieldMixin!();
}

class ExecuterContext:CalcContext{
    TemplateExecuter input;
    TemplateHandler templateH;

    this(TemplateExecuter input,CalculationInstance cInstance,char[] className,char[] contextId){
        super(cInstance,contextId);
        this.input=input;
        templateH=new TemplateHandler(input.templateDirectory(),cInstance.baseDirectory());
        input.addFullSubs(templateH.subs);
        
        templateH.subs["templateDirectory"]=input.templateDirectory.toString;
        templateH.subs["workingDirectory"]=templateH.targetDir.toString;
        
        templateH.evalTemplates(0,true);
        cInstance.execCmd(input.setupCommand());
        sout("pippo_a6..\n");
        if (input.startConfig is null || input.startConfig.config is null){
            throw new Exception("Error: startConfiguration in field "~input.myFieldName~" should be set to a valid configuration",__FILE__,__LINE__);
        }
        auto pSys=input.startConfig.config.particleSysReal();
        
        sout("pippo_a7\n");
        if (pSys.nCenter is null)
            pSys.nCenter=nCenter;
        else
            _nCenter=pSys.nCenter;
        sout("pippo_a8\n");
        
        templateH.longSubs["coord.xyz"]=&writeXyz;
        templateH.longSubs["turboCoord"]=&writeTurboCoord;
        
        sout("pippo_a9\n");
        templateH.evalTemplates(0,true);
        sout("pippo_a10\n");
    }
    
    void writeXyz(CharSink s){
        mixin(withPSys("WriteOut.writeXyz(s,pSys.sysStruct,pSys.dynVars.pos,pSys.name);"));
    }

    void writeTurboCoord(CharSink s){
        mixin(withPSys("WriteOut.writeTurboCoord(s,pSys.sysStruct,pSys.dynVars.pos);"));
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        templateH.evalTemplates(changeLevel,input.overwriteUnchangedPaths);
        cInstance.execCmd(input.commandFor(updateE,updateF,changeLevel));
        collectEF(updateE,updateF);
        maxChange=0;
        changeLevel=ChangeLevel.SmoothPosChange;
    }
    
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        assert(0,"to implement in subclasses");
    }
    
    void activate(){
        cInstance.activatedContext(this);
    }
    void deactivate(){
        cInstance.deactivatedContext(this);
    }
    void destroy(){
        cInstance.destroyedContext(this);
    }
}


class TemplateHandler{
    char[][char[]] subs;
    OutWriter[char[]] longSubs;
    VfsFolder sourceDir;
    VfsFolder targetDir;
    size_t maxKeywordLen=128;
    bool shouldWriteReplacementsDict;
    char[] rDictFilename;
    
    this(VfsFolder sourceDir,VfsFolder targetDir,char[][char[]] subs=null){
        this.sourceDir=sourceDir;
        this.targetDir=targetDir;
        this.subs=subs;
    }
    /// writes the current replacements dictionary to the given sink
    void writeReplacementsDict(CharSink outF){
        void writeEscaped(char[] str){
            size_t wrote=0;
            foreach(i,c;str){
                switch (c){
                case '\\','\"':
                    outF(str[wrote..i]);
                    outF("\\");
                    wrote=i;
                    break;
                default:
                    break;
                }
            }
            outF(str[wrote..$]);
        }
        outF("{");
        bool cont=false;
        foreach(k,v;subs){
            if (cont) outF(",");
            outF("\n  \"");
            cont=true;
            writeEscaped(k);
            outF("\":\"");
            writeEscaped(v);
            outF("\"");
        }
        foreach(k,v;longSubs){
            if (cont) outF(",");
            outF("\n  \"");
            cont=true;
            writeEscaped(k);
            outF("\":\"");
            v(&writeEscaped);
            outF("\"");
        }
        outF("\n}\n");
    }
    /// evaluates the templates
    void evalTemplates(int level,bool overwriteUnchangedPaths=false){
        if (shouldWriteReplacementsDict){
            auto rDictF=targetDir.file(rDictFilename).create();
            auto outF=strDumper(rDictF.output);
            writeReplacementsDict(outF);
            rDictF.output.close();
        }
        foreach (f;sourceDir.tree.catalog()){
            if (f.name.length>7 && f.name[$-7..$-1]==".templ" &&
                f.name[$-1]>='0' && f.name[$-1]<='9')
            {
                auto lFile=cast(int)(f.name[$-1]-'0');
                auto newName=f.name()[0..$-7]; // should I use toString??
                if ((lFile==level || (lFile<level && 
                     (!targetDir.file(f.toString()[0..$-8]~(cast(char)('0'+level))).exists)) &&
                    !targetDir.file(newName).exists()))
                {
                    auto fIn=f.input;
                    scope(exit){
                        fIn.close();
                    }
                    auto newF=targetDir.file(newName).create();
                    auto fOut=newF.output;
                    scope(exit){
                        fOut.close();
                    }
                    makeSubs(strReaderHandler(fIn),strDumper(fOut));
                }
                // 0: all changed, 1: only pos changed, 2: small pos change, 3: extrapolable pos change
            } else if (overwriteUnchangedPaths || level==0 ||
                (!targetDir.file(f.toString()).exists()))
            {
                auto newF=targetDir.file(f.name).create(f.input);
                auto output=newF.output;
            }
        }
    }
    
    /// possibly replaces name with a template expansion into outF
    bool maybeReplace(char[] name,CharSink outF){
        auto sub=name in subs;
        if (sub !is null){
            assert((name in longSubs)is null,"double entry");
            outF(*sub);
            return true;
        }
        auto sub2=name in longSubs;
        if (sub2 !is null){
            (*sub2)(outF);
            return true;
        }
        return false;
    }
    
    /// writes the content of inF to outF, performing substitutions on it
    void makeSubs(bool delegate(CharReader) inF,CharSink outF){
        inF(delegate size_t(char[]data, SliceExtent slice,out bool iterate){
            iterate=true;
            for(size_t i=0;i<data.length;++i){
                if (data[i]=='['){
                    if (i==0){
                        for(size_t j=i;j<data.length;++j){
                            if(data[j]=='['){
                                outF(data[0..j]);
                                return j;
                            }
                            if(data[j]==']'){
                                if (!maybeReplace(data[i..(j+1)],outF)) outF(data[i..(j+1)]);
                                return j+1;
                            }
                        }
                        if(data.length<maxKeywordLen){
                            switch (slice){
                            case SliceExtent.Partial:
                                return Eof;
                            case SliceExtent.Maximal:
                                outF(data);
                                throw new BIOException("unexpected SliceExtent.Maximal with less than maxKeywordLen size",__FILE__,__LINE__);
                            case SliceExtent.ToEnd:
                                iterate=false;
                                outF(data);
                                return data.length;
                            default: assert(0);
                            }
                        }
                        return data.length;
                    } else {
                        outF(data[0..i]);
                        return i;
                    }
                }
            }
            switch (slice){
            case SliceExtent.Partial:
                if (data.length==0)
                    return Eof;
                break;
            case SliceExtent.Maximal:
                if (data.length==0)
                    throw new BIOException("unexpected SliceExtent.Maximal with 0 size",__FILE__,__LINE__);
                break;
            case SliceExtent.ToEnd:
                iterate=false;
                return data.length;
            default: assert(0);
            }
            outF(data);
            return data.length;
        });
    }
}