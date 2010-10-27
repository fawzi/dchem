module dchem.calculator.FileCalculator;
import blip.io.BasicIO;
import blip.io.StreamConverters;
import tango.io.vfs.FileFolder;
import blip.time.Time;
import blip.time.Clock;
import WriteOut=dchem.input.WriteOut;
import dchem.calculator.Calculator;
import dchem.input.RootInput;
import blip.sync.UniqueNumber;
import blip.serialization.Serialization;
import tango.sys.Process;
import blip.container.GrowableArray;
import Path=tango.io.Path;
import blip.io.Console;
import dchem.util.ExecCmd;
import dchem.calculator.ProcContext;
import dchem.input.ReadCoordFile;

class TemplateExecuter: Method {
    InputField superTemplate;
    InputField startConfig;
    char[] templateDir;
    char[] setupCmd;
    char[] stopCmd;
    char[][char[]] subs;
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
        auto s=dumper(sink);
        if (superTemplate !is null){
            if ((cast(TemplateExecuter)superTemplate.content)is null){
                s("Error: superTemplate should be of type dchem.TemplateExecuter in field "~myFieldName);
                res=false;
            }
        }
        return res;
    }
    
    CalculationContext getCalculator(bool wait,ubyte[]history){
        assert(0,"to be implemented by subclasses");
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
        return "";
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
    startConfig: the initial configuration
    templateDir: where to find the definition of the template (tipically a directory)
    subs: keyword and their substitutions to apply to the templates (as dictionary string -> string)
    setupCmd: a command that should be executed to set up the context
    stopCmd: a command that should be executed to stop the context
    writeReplacementsDict: if a dictionary with the replacements performed should be written (default is false)
    overwriteUnchangedPaths: if paths that are already there should be overwitten (default is false)
    maxContexts: maximum number of contexts per calculation instance (default is 1)
    `));
    mixin printOut!();
    mixin myFieldMixin!();
    
    /// drops the history associated with the given key
    void dropHistory(ubyte[]history){}
    /// clears all history
    void clearHistory(){}
}

class CmdTemplateExecuter:TemplateExecuter {
    char[] executeInitialE, executeStructChangeE, executePosChangeE, executeSmallPosChangeE, executeSmoothPosChangeE, executeDefaultE;
    char[] executeInitialEF, executeStructChangeEF, executePosChangeEF, executeSmallPosChangeEF, executeSmoothPosChangeEF, executeDefaultEF;
    
    bool verify(CharSink sink){
        bool res=super.verify(sink);
        return res;
    }
    
    CalculationContext getCalculator(bool wait,ubyte[]history){
        assert(0,"to be implemented by subclasses");
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
    mixin(serializeSome("dchem.CmdTemplateExecuter",`
    executeInitialE: command executed to calculate the first energy, set to NONE to deactivate
    executeStructChangeE: command executed to calculate the energy when the structure is changed, set to NONE to deactivate
    executePosChangeE: command executed to calculate the energy when only positions changed, set to NONE to deactivate
    executeSmallPosChangeE: command executed to calculate the energy when positions changed by a small amount, set to NONE to deactivate
    executeSmoothPosChangeE: command executed to calculate the energy when positions changed smoothly, set to NONE to deactivate
    executeDefaultE: default command for all execute*E commands (if they are empty)
    executeInitialEF: command executed to calculate the first energy and forces, set to NONE to deactivate
    executeStructChangeEF: command executed to calculate the energy and forces when the structure is changed, set to NONE to deactivate
    executePosChangeEF: command executed to calculate the energy and forces when only positions changed, set to NONE to deactivate
    executeSmallPosChangeEF: command executed to calculate the energy and forces when positions changed by a small amount, set to NONE to deactivate
    executeSmoothPosChangeEF: command executed to calculate the energy and forces when positions changed smoothly, set to NONE to deactivate
    executeDefaultEF: default command for all execute*EF commands (if they are empty)
    `));
}

class ExecuterContext:CalcContext{
    TemplateExecuter input;
    TemplateHandler templateH;
    Process opInProgress;

    VfsFolder baseDir(){
        return templateH.targetDir;
    }
    Method method(){
        return input;
    }
    void execCmd(char[] cmd,CharSink log=sout.call){
        if (cmd.length>0 && cmd!="NONE"){
            char[256] buf;
            auto arr=lGrowableArray(buf,0);
            dumper(&arr)("executing ")(cmd)("\n");
            log(arr.data);
            arr.clearData();
            opInProgress=new Process(cmd);
            auto bDir=baseDir.toString();
            if (bDir.length>0){
                opInProgress.workDir=bDir;
            }
            int status;
            log(getOutput(opInProgress,status));
            if (status!=0){
                arr("command failed with status ");
                writeOut(&arr.appendArr,status);
                throw new Exception(arr.data.dup,__FILE__,__LINE__);
            }
        }
    }

    this(TemplateExecuter input,char[] contextId){
        super(contextId);
        this.input=input;
        templateH=new TemplateHandler(input.templateDirectory(),new FileFolder(ProcContext.instance.baseDirectory.toString()~"/"~contextId,true)); // to fix
        input.addFullSubs(templateH.subs);
        
        templateH.subs["templateDirectory"]=input.templateDirectory.toString;
        templateH.subs["workingDirectory"]=templateH.targetDir.toString;
        
        templateH.evalTemplates(0,true);
        execCmd(input.setupCommand());
        if (input.startConfig is null || cast(Config)input.startConfig.contentObj is null){
            throw new Exception("Error: startConfiguration in field "~input.myFieldName~" should be set to a valid configuration",__FILE__,__LINE__);
        }
        auto pSys=(cast(Config)input.startConfig.contentObj).particleSysReal();
        
        if (pSys.nCenter is null)
            pSys.nCenter=nCenter;
        else
            _nCenter=pSys.nCenter;
        _pSysReal=pSys;
        
        templateH.longSubs["coord.xyz"]=&writeXyz;
        templateH.longSubs["turboCoord"]=&writeTurboCoord;
        
        templateH.evalTemplates(0,true);
    }
    
    void writeXyz(CharSink s){
        mixin(withPSys("WriteOut.writeXyz(s,pSys.sysStruct,pSys.dynVars.x.pos,pSys.name);"));
    }

    void writeTurboCoord(CharSink s){
        mixin(withPSys("WriteOut.writeTurboCoord(s,pSys.sysStruct,pSys.dynVars.x.pos);"));
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        templateH.evalTemplates(changeLevel,input.overwriteUnchangedPaths);
        execCmd(input.commandFor(updateE,updateF,changeLevel));
        collectEF(updateE,updateF);
        maxChange=0;
        changeLevelSet=ChangeLevel.SmoothPosChange;
    }
    
    override void stop(){
        volatile auto op=opInProgress;
        if (op!is null && op.isRunning)
            op.kill();
    }
    
    /// should collect the newly calculated energy
    void collectEF(bool updateE=true,bool updateF=true){
        assert(0,"to implement in subclasses");
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