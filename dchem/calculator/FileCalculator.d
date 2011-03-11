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
import tango.io.stream.Data;
import blip.container.GrowableArray;
import Path=tango.io.Path;
import blip.io.Console;
import dchem.util.ExecCmd;
import dchem.calculator.ProcContext;
import dchem.input.ReadCoordFile;
import blip.parallel.mpi.MpiModels;
import dchem.Common;
import blip.io.FileStream;
import blip.core.Array;

/// helper for defining things to log to special files
struct EvalLog{
    string targetFile;
    string format;
    mixin(serializeSome("dchem.EvalLog",`
    targetFile: file where to log the things
    format: the format (xyz,xyzForces,energy,jsonDynMddx,sbinDynMddx,jsonDynX,sbinDynX,jsonFull,sbinFull)`));
    mixin printOut!();
}

/// an executer that uses a template directory, and is command based
class TemplateExecuter: Method {
    InputField superTemplate;
    InputField startConfig;
    char[] templateDir;
    char[] setupCmd;
    char[] setupCtxCmd;
    char[] stopCtxCmd;
    char[][char[]] subs;
    bool writeReplacementsDict;
    bool overwriteUnchangedPaths;
    bool makeReplacementsInCommands=true;
    bool ignoreSetupExitStatus;
    bool ignoreSetupCtxExitStatus;
    EvalLog[] onELog=[{targetFile:"log.energies",format:"energy"},
        {targetFile:"log.pos",format:"xyz"}];
    EvalLog[] onFLog=[{targetFile:"log.forces",format:"xyzForces"}];
    VfsFolder _templateDirectory;
    CharSink logger;
    
    this(){
    }
    
    TemplateExecuter sTemplate(){
        if (superTemplate!is null) {
            return superTemplate.contentT!(TemplateExecuter)();
        }
        return null;
    }
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
            if ((cast(TemplateExecuter)(superTemplate.contentObj()))is null){
                s("Error: superTemplate should be of type dchem.TemplateExecuter in field "~myFieldName);
                res=false;
            }
        }
        foreach(l;onELog){
            if (l.targetFile.length==0){
                dumper(sink)("missing targetFile in onELog in field ")(myFieldName)("\n");
                res=false;
            }
            if (find(WriteOut.writeConfigFormats,l.format)==WriteOut.writeConfigFormats.length){
                auto w=dumper(sink);
                w("format in onELog has to be one of ");
                foreach(i,f;WriteOut.writeConfigFormats){
                    if (i!=0) w(", ");
                    w(f);
                }
                w(" and not '")(l.format)("' in field ")(myFieldName)("\n");
                res=false;
            }
        }
        foreach(l;onFLog){
            if (l.targetFile.length==0){
                dumper(sink)("missing targetFile in onFLog in field ")(myFieldName)("\n");
                res=false;
            }
            if (find(WriteOut.writeConfigFormats,l.format)==WriteOut.writeConfigFormats.length){
                auto w=dumper(sink);
                w("format in onFLog has to be one of ");
                foreach(i,f;WriteOut.writeConfigFormats){
                    if (i!=0) w(", ");
                    w(f);
                }
                w(" and not '")(l.format)("' in field ")(myFieldName)("\n");
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
            sTemplate().setupCommand();
        }
        return setupCmd;
    }
    char[] setupCtxCommand(){
        if (setupCtxCmd.length==0 && superTemplate!is null){
            sTemplate().setupCtxCommand();
        }
        return setupCtxCmd;
    }
    char[] stopCtxCommand(){
        if (stopCtxCmd.length==0 && superTemplate!is null){
            sTemplate().stopCtxCommand();
        }
        return stopCtxCmd;
    }
    char[] commandFor(bool energy,bool force,ChangeLevel changeLevel,Real change){
        return "";
    }
    void addFullSubs(ref char[][char[]] sub){
        assert(sub!is subs||sub is null,"cannot add directly to own subs");
        if (superTemplate!is null){
            sTemplate().addFullSubs(sub);
        }
        foreach (k,v;subs){
            sub[k]=v;
        }
    }
    void setup(LinearComm pEnv,CharSink log){
        if (sTemplate()!is null){
            sTemplate().setup(pEnv,log);
        }
        logger=log;
        if (setupCommand.length>0 && setupCommand!="NONE"){
            auto templateH=new TemplateHandler(templateDirectory(),new FileFolder(ProcContext.instance.baseDirectory.toString(),true)); // to fix
            addFullSubs(templateH.subs);
            templateH.shouldWriteReplacementsDict=writeReplacementsDict;
            templateH.subs["templateDirectory"]=templateDirectory.toString;
            templateH.subs["workingDirectory"]=templateH.targetDir.toString;
            auto opInProgress=templateH.processForCmd(setupCommand,log);
            templateH.exec(opInProgress,log,ignoreSetupExitStatus);
        }
    }
    // serialization stuff
    mixin(serializeSome("dchem.TemplateExecuter",
    `superTemplate: a template where to take default values
    startConfig: the initial configuration
    templateDir: where to find the definition of the template (tipically a directory)
    subs: keyword and their substitutions to apply to the templates (as dictionary string -> string)
    setupCmd: a command that should be executed when the method is setup
    setupCtxCmd: a command that should be executed to set up the context
    stopCtxCmd: a command that should be executed to stop the context
    writeReplacementsDict: if a dictionary with the replacements performed should be written (default is false)
    overwriteUnchangedPaths: if paths that are already there should be overwitten (default is false)
    makeReplacementsInCommands: if replacements should be performed on the commands to be executed (default is true)
    ignoreSetupExitStatus: if the exit status of the setup command should be ignored (false)
    ignoreSetupCtxExitStatus: if the exit status of the context setup command should be ignored (false)
    onELog: what to log for each energy evaluation (by default energy and positions ad xyz)
    onFLog: what to log for each force evaluation (by default xyzForces)
    `));
    mixin printOut!();
    mixin myFieldMixin!();
    
    /// drops the history associated with the given key
    void dropHistory(ubyte[]history){}
    /// clears all history
    void clearHistory(){}
}

class CmdTemplateExecuter:TemplateExecuter {
    char[] executeE;
    char[] executeEF;
    char[] executeF0;
    
    bool verify(CharSink sink){
        bool res=super.verify(sink);
        return res;
    }
    CalculationContext getCalculator(bool wait,ubyte[]history){
        assert(0,"to be implemented by subclasses");
    }
    char[] commandFor(bool energy,bool force,ChangeLevel changeLevel,Real change){
        char[] res;
        if (force==false){
            if (energy==false) return "NONE";
            res=executeE;
        } else if (energy==false && changeLevel==ChangeLevel.NoChange){
            res=executeF0;
        }
        if (res.length==0) res=executeEF;
        if (res.length==0 && superTemplate !is null){
            sTemplate().commandFor(energy,force,changeLevel,change);
        }
        return res;
    }
    void addFullSubs(ref char[][char[]] sub){
        assert(sub!is subs||sub is null,"cannot add directly to own subs");
        if (superTemplate!is null){
            sTemplate().addFullSubs(sub);
        }
        foreach (k,v;subs){
            sub[k]=v;
        }
    }
    // serialization stuff
    mixin(serializeSome("dchem.CmdTemplateExecuter",`
    executeE: command executed to calculate the energy alone (optional)
    executeEF: command executed to calculate the energy and forces (this is required)
    executeF0: command executed do calculate the forces immediately after the energy (optional)`));
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
    void execCmd(char[] cmd,CharSink log=null,bool ignoreExitStatus=false){
        if (log is null) log=logger;
        if (cmd.length>0 && cmd!="NONE"){
            opInProgress=templateH.processForCmd(cmd,log);
            if (opInProgress is null) throw new Exception("could not create process for command "~cmd,__FILE__,__LINE__);
            templateH.exec(opInProgress,log,ignoreExitStatus);
        }
    }
    
    static TemplateHandler initialTH(TemplateExecuter input,char[] contextId){
        auto templateH=new TemplateHandler(input.templateDirectory(),new FileFolder(ProcContext.instance.baseDirectory.toString()~"/"~contextId,true)); // to fix
        input.addFullSubs(templateH.subs);
        templateH.shouldWriteReplacementsDict=input.writeReplacementsDict;
        templateH.subs["templateDirectory"]=input.templateDirectory.toString;
        templateH.subs["workingDirectory"]=templateH.targetDir.toString;
        return templateH;
    }

    this(TemplateExecuter input,char[] contextId,TemplateHandler th=null){
        super(contextId,input.logger);
        this.input=input;
        templateH=th;
        if (th is null){
            templateH=initialTH(input,contextId);
        }
        templateH.evalTemplates(0,true);
        execCmd(input.setupCtxCommand(),null,input.ignoreSetupCtxExitStatus);
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
    
    override void desc(CharSink s){
        dumper(s)("Executer '")(contextId())("' generated by field '")(input.myFieldName)("'");
    }
    
    void logAfterEval(bool updatedE,bool updatedF){
        if (updatedE){
            foreach(l;input.onELog){
                auto f=outfileBin(Path.join(templateH.targetDir.toString,l.targetFile),WriteMode.WriteAppend);
                scope(exit){ f.flush(); f.close(); }
                mixin(withPSys(`
                WriteOut.writeConfig(f,pSys,l.format,externalRef);
                `));
            }
        }
        if (updatedF){
            foreach(l;input.onFLog){
                auto f=outfileBin(Path.join(templateH.targetDir.toString,l.targetFile),WriteMode.WriteAppend);
                scope(exit){ f.flush(); f.close(); }
                mixin(withPSys(`
                WriteOut.writeConfig(f,pSys,l.format,externalRef);
                `));
            }
        }
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        templateH.subs["evalE"]=(updateE?"T":"F");
        templateH.subs["evalF"]=(updateF?"T":"F");
        templateH.evalTemplates(changeLevel,input.overwriteUnchangedPaths);
        execCmd(input.commandFor(updateE,updateF,changeLevel,maxChange));
        try{
            collectEF(updateE,updateF);
        } catch (Exception e){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("Error trying to collect energy/forces in ")(&this.desc);
            }),__FILE__,__LINE__,e);
        }
        templateH.subs["evalE"]="-";
        templateH.subs["evalF"]="-";
        maxChange=0;
        changeLevelSet=ChangeLevel.NoChange;
        logAfterEval(updateE,updateF);
    }
    
    override void stop(){
        volatile auto op=opInProgress;
        auto stopC=input.stopCtxCmd;
        if (stopC!is null && op !is null){
            execCmd(stopC,null,true);
        } else if (op!is null && op.isRunning){
            op.kill();
        }
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
    char[] rDictFilename="repl.dict";
    
    this(VfsFolder sourceDir,VfsFolder targetDir,char[][char[]] subs=null){
        this.sourceDir=sourceDir;
        this.targetDir=targetDir;
        this.subs=subs;
    }
    /// returns a Process that executes the given command
    Process processForCmd(char[] cmd,CharSink log=sout.call,bool makeReplacementsInCommands=true){
        Process opInProgress;
        if (cmd.length>0 && cmd!="NONE"){
            char[256] buf;
            auto arr=lGrowableArray(buf,0);
            dumper(&arr)("will execute `")(cmd)("`\n");
            log(arr.data);
            arr.clearData();
            auto cmdP=cmd;
            if (cmd.length!=0 && makeReplacementsInCommands){
                makeSubs(delegate bool(CharReader r){
                    bool res=false;
                    while(cmdP.length!=0){
                        bool iterate=false;
                        auto l=r(cmdP,SliceExtent.ToEnd,iterate);
                        if (l!=0){
                            if (l>cmdP.length) throw new Exception("requested more than available",__FILE__,__LINE__);
                            cmdP=cmdP[l..$];
                            res=true;
                        }
                        if (!iterate) {
                            return res;
                        }
                    }
                    return res;
                },&arr.appendArr);
                opInProgress=new Process(arr.data.dup);
                arr(" after substitutions\n");
                log(arr.data);
            } else {
                opInProgress=new Process(cmdP);
            }
            auto bDir=targetDir.toString();
            if (bDir.length>0){
                opInProgress.workDir=bDir;
            }
            opInProgress.copyEnv=true;
        }
        return opInProgress;
    }
    /// executes the given Process
    void exec(Process opInProgress,CharSink log=sout.call,bool ignoreExitStatus=false){
        int status;
        log(getOutput(opInProgress,status));
        if (status!=0){
            char[256] buf;
            auto arr=lGrowableArray(buf,0,GASharing.Local);
            arr("command failed with status ");
            writeOut(&arr.appendArr,status);
            if (ignoreExitStatus){
                arr("\n");
                log(arr.data);
            } else {
                throw new Exception(arr.takeData,__FILE__,__LINE__);
            }
        }
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
        subs["changeLevel"]=collectAppender(delegate void(CharSink s){
            writeOut(s,level);
        });
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
                if (lFile==level || (lFile>level && 
                        (!targetDir.file(f.toString()[0..$-8]~(cast(char)('0'+level))).exists)) ||
                    (!targetDir.file(newName).exists()))
                {
                    auto fIn=new DataInput(f.input);
                    scope(exit){
                        fIn.close();
                    }
                    auto newF=targetDir.file(newName).create();
                    auto fOut=newF.output;
                    scope(exit){
                        fOut.flush();
                        fOut.close();
                    }
                    makeSubs(strReaderHandler(fIn),strDumper(fOut));
                }
                // 0: all changed, 1: only pos changed, 2: small pos change, 3: extrapolable pos change
            } else if (overwriteUnchangedPaths || level==0 ||
                (!targetDir.file(f.toString()).exists()))
            {
		auto fIn=f.input;
		scope(exit){
		    fIn.close();
		}
		auto newF=targetDir.file(f.name).create(fIn);
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
                        for(size_t j=i+1;j<data.length;++j){
                            if(data[j]=='['){
                                outF(data[0..j]);
                                return j;
                            }
                            if(data[j]==']'){
                                if (!maybeReplace(data[(i+1)..j],outF)) outF(data[i..(j+1)]);
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
            default: assert(0);
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
            }
            outF(data);
            return data.length;
        });
    }
}
