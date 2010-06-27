/// describes the root input file
module dchem.input.RootInput;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import blip.serialization.Serialization;
import blip.text.TextParser;
import blip.io.BasicIO;
import blip.io.Console;
import dchem.input.ReadIn;
import dchem.Common;
import blip.util.TangoLog;
import blip.core.Variant;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import tango.io.FilePath;
import dchem.input.ReadIn2PSys;
import blip.util.NotificationCenter;
import dchem.sys.PIndexes;
import blip.BasicModels;
import dchem.sys.Constraints;
import tango.io.vfs.model.Vfs;
//import tango.sys.Process;

/// represents a task that can be sent over to another computer
interface RemoteTask:Serializable{
    void execute(Variant args);
    void stop();
}

/// an element of the input
interface InputElement:Serializable{
    /// returns the field of this object
    InputField myField();
    /// called by the containing InputField after deserialization
    void myField(InputField);
    /// should check the input (called after full deserialization)
    bool verify(CharSink logger);
}

enum InstanceGetFlags{
    ReuseCache=1,    /// reuse instances that are in the cache
    NoAlloc=2,       /// don't alloc a new instance
    Wait=4,          /// wait until an instance is available (when one hits the maximum number of instances)
    ReuseWait=5      /// ReuseCache|Wait
}

/// calculator setup (chooses the method to perform each calculation)
interface Method:InputElement{
    /// gets a calculator to perform calculations with this method, if possible reusing the given history
    /// if wait is true waits until a context is available
    CalculationContext getCalculator(bool wait,ubyte[]history);
    /// drops the history with the given id
    void dropHistory(ubyte[]history);
    /// clears all history
    void clearHistory();
    /// url to access this from other processes
    char[] exportedUrl();
}


/// configuration
interface Config:InputElement{
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem();
    /// the configuration as ParticleSys (might be cached)
    ParticleSys!(Real) particleSysReal();
    /// the configuration as low precision ParticleSys (might be cached)
    ParticleSys!(LowP) particleSysLowP();
    /// drops the cached readSystem/particleSys
    void clear();
}

/// amount of change since the last calculation in the context
enum ChangeLevel{
    FirstTime=0,
    AllChanged=1,
    PosChanged=2,
    SmallPosChange=3,
    SmoothPosChange=4
}
/// represent a calculation that might have been aready partially setup, in particular the
/// number of elements,... cannot change
interface CalculationContext{
    /// unique identifier for this context
    char[] contextId();
    /// the particle system based on real numbers (will be null if pSysLowP is valid)
    ParticleSys!(Real) pSysReal();
    /// the particle system based on low precision numbers (will be null if pSysReal is valid)
    ParticleSys!(LowP) pSysLowP();
    /// the system struct
    SysStruct sysStruct();
    /// notification central of the current particle system & context.
    /// will receive activation, deactivation & destruction notifications
    NotificationCenter nCenter();
    /// history of the previous positions
    HistoryManager!(LowP) posHistory();
    /// change level since the last calculation
    ChangeLevel changeLevel();
    /// sets the change level
    void changeLevel(ChangeLevel);
    /// decreases the change level to at least changeLevel
    void changedDynVars(ChangeLevel changeLevel,Real diff);
    /// the total potential energy
    Real potentialEnergy();
    /// changes the position (utility method)
    void pos(SegmentedArray!(Vector!(Real,3)) newPos);
    /// gets the positions (utility method)
    SegmentedArray!(Vector!(Real,3)) pos();
    /// changes the velocities (utility method)
    void dpos(SegmentedArray!(Vector!(Real,3)) newDpos);
    /// gets the velocities (utility method)
    SegmentedArray!(Vector!(Real,3)) dpos();
    /// changes the forces (utility method)
    void mddpos(SegmentedArray!(Vector!(Real,3)) newDdpos);
    /// gets the forces (utility method)
    SegmentedArray!(Vector!(Real,3)) mddpos();
    /// updates energy and/or forces
    void updateEF(bool updateE=true,bool updateF=true);
    /// called automatically after creation, but before any energy evaluation
    /// should be called before working again with a deactivated calculator
    void activate();
    /// call this to possibly compress the context and empty caches (i.e. before a long pause in the calculation)
    void deactivate();
    /// call this to gives back the context (after all calculations with this are finished) might delete or reuse it
    void giveBack();
    /// tries to stop a calculation in progress, recovery after this is not possible (only giveBack can be called)
    void stop();
    /// a place where to put all active constraints
    MultiConstraint constraints();
    /// method that has generated this context
    Method method();
    /// stores the history somewhere and returns an id to possibly recover that history at a later point
    /// this is just an optimization, it does not have to do anything. If implemented then
    /// method.getCalculator, .dropHistory and .clearHistory have to be implemented accordingly
    ubyte[]storeHistory();
    /// url to access this from other processes
    char[] exportedUrl();
}

/// templatized way to extract the pSys from a CalculationContext
ParticleSys!(T) pSysT(T)(CalculationContext c){
    static if (is(Real==T)){
        return c.pSysReal;
    } else static if(is(LowP==T)){
        return c.pSysLowP;
    } else {
        static assert("requested a pSys "~T.stringof~" which is neither Real ("~
            Real.stringof~") nor LowP ("~LowP.stringof~")");
    }
}

/// a calculation method
interface Sampler:InputElement{
    void run();
    void stop();
}

/// template to mixin to give myField & co
template myFieldMixin(){
    InputField _myField;
    InputField myField() { return _myField; }
    void myField(InputField i){ _myField=i; }
    char[] myFieldName(){
        if (_myField is null){
            return "*null*";
        } else {
            return _myField.name;
        }
    }
}

/// an root input entry
/// serializes as reference, unserializes normally
class InputField:InputElement{
    enum TypeId{
        Configuration,
        Method,
        Executer,
        Sampler,
        Reference,
        InputField,
    }
    char[] name;
    InputElement content;
    TypeId typeId;
    mixin myFieldMixin!();
    
    /// constructor just for unserialization
    this(){
        this("dummy");
    }
    /// constructor
    this(char[]name,TypeId typeId=TypeId.Reference,InputElement content=null){
        this.name=name;
        this.typeId=typeId;
        this.content=content;
    }

    static char[]typeStr(TypeId t){
        switch (t){
            case TypeId.Configuration: return "config";
            case TypeId.Method: return "method";
            case TypeId.Executer: return "executer";
            case TypeId.Sampler: return "sampler";
            case TypeId.Reference: return "ref";
            case TypeId.InputField: return "inputField";
            default: assert(0,"invalid typeId");
        }
    }
    static TypeId strType(char[] val){
        switch (val){
            case "config":  return TypeId.Configuration;
            break;
            case "method":  return TypeId.Method;
            break;
            case "sampler": return TypeId.Sampler;
            break;
            case "ref":     return TypeId.Reference;
            break;
            case "inputField": return TypeId.InputField;
            break;
            default: assert(0,"invalid typeId");
        }
    }
    /// sets typeId depending on content
    static TypeId objType (Serializable o){
        if (o is null) return TypeId.Reference;
        auto oo=cast(Object)o;
        if (cast(Config)oo) return TypeId.Configuration;
        if (cast(Method)oo) return TypeId.Method;
        if (cast(Sampler)oo) return TypeId.Sampler;
        if (cast(InputField)oo) return TypeId.InputField;
        throw new Exception("could not find type of object",__FILE__,__LINE__);
    }
    Config config(){
        return cast(Config)cast(Object)content;
    }
    Method method(){
        return cast(Method)cast(Object)content;
    }
    Sampler sampler(){
        return cast(Sampler)cast(Object)content;
    }
    InputField inputField(){
        return cast(InputField)cast(Object)content;
    }
    // serialization stuff
    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(this))("Ref");
        metaI.addFieldOfType!(char[])("name","name of the reference");
    }
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    void preSerialize(Serializer s){ }
    void postSerialize(Serializer s){ }
    void serialize(Serializer s){
        s.field(metaI[0],name);
        auto rInp0="RootInput" in s.context;
        RootInput rInp;
        if (rInp0 is null){
            Log.lookup("blip.parallel.smp").warn("serializing InputField without a RootInput in serializer context",
                __FILE__,__LINE__);
            rInp=new RootInput();
            s.context["RootInput"]=rInp;
        } else {
            rInp=rInp0.get!(RootInput)();
        }
        rInp.addName(this);
    }
    Serializable preUnserialize(Unserializer s){
        return this;
    }
    /// unserialize an object
    void unserialize(Unserializer s){
        s.field(metaI[0],name);
        auto rInp0="RootInput" in s.context;
        RootInput rInp;
        if (rInp0 is null){
            Log.lookup("blip.parallel.smp").warn("warning: unserializing InputField without a RootInput in serializer context",
                __FILE__,__LINE__);
            rInp=new RootInput();
            s.context["RootInput"]=rInp;
        } else {
            rInp=rInp0.get!(RootInput)();
        }
    }
    Serializable postUnserialize(Unserializer s){
        return s.context["RootInput"].get!(RootInput)().makeUnique(this);
    }
    /// unify with another value
    void unify(InputField o,bool overwrite=false){
        if (o!is null && o.typeId!=TypeId.Reference){
            if (typeId==TypeId.Reference){
                typeId=o.typeId;
                content=o.content;
            } else {
                if (o.typeId!=typeId){
                    throw new Exception("unification of InputFields of different types, duplicate entry for "~name~"?",__FILE__,__LINE__);
                }
                if (content!is o.content && !overwrite){
                    throw new Exception("unification of InputFields changed value and overwrite is false",__FILE__,__LINE__);
                }
                content=o.content;
            }
        }
    }
    bool verify(CharSink log){
        return true;
    }
}

/// represents the main input of the program
class RootInput{
    InputField[char[]] knownNames;
    char[][] newKnownNames;
    bool rewriteNamesOk;
    
    /// reads the input from the given TextParser
    bool readInput(TextParser!(char) t,CharSink errorLog){
        scope js=new JsonUnserializer!(char)(t);
        js.sloppyCommas=true;
        js.context["RootInput"]=Variant(this);
        t.newlineIsSpace=true;
        t.skipWhitespace();
        if (t.getSeparator()!="{"){
            t.parseError("unexpected start expecting '{'",__FILE__,__LINE__);
        }
        while (true){
            t.newlineIsSpace=true;
            t.skipWhitespace();
            auto sep=t.getSeparator();
            if (sep.length!=0){
                if (sep=="}") break;
                if (sep!=",") {
                    t.parseError("uexpected separator '"~sep~"'",__FILE__,__LINE__);
                }
            }
            if (!t.next(&t.scanString)){
                sep=t.getSeparator();
                if (sep=="}"){
                    break;
                } else if (sep.length==0) {
                    t.parseError("unexpected EOF, expecting '}'",__FILE__,__LINE__);
                } else {
                    t.parseError("unexpected separator '"~sep~"'",__FILE__,__LINE__);
                }
            }
            char[] name=maybeUnescape(t.get.dup);
            if (t.getSeparator()!=":"){
                throw new Exception("unexpected separator",__FILE__,__LINE__);
            }
            // check if it is already there
            auto inputF=new InputField(name);
            inputF=makeUnique(inputF);
            InputElement content;
            js(content);
            if (inputF.content!is null && !rewriteNamesOk){
                throw new Exception("duplicated name '"~name~"' in input",__FILE__,__LINE__);
            }
            inputF.content=content;
            inputF.content.myField(inputF);
            inputF.typeId=InputField.objType(inputF.content);
        }
        bool inputOk=true;
        foreach (k,v;knownNames){
            if (v.content is null){
                errorLog("Error: field "~k~" is null\n");
                inputOk=false;
            } else {
                inputOk=inputOk && v.content.verify(errorLog);
            }
        }
        return inputOk;
    }
    
    /// reads the input from the given TextParser
    void writeInput(CharSink s){
        auto nNames=knownNames.keys;
        newKnownNames=null;
        scope js=new JsonSerializer!(char)(s);
        js.context["RootInput"]=Variant(this);
        while (true){
            foreach(k;nNames) {
                js.handlers(k);
                s(":");
                js(knownNames[k].content);
            }
            nNames=newKnownNames;
            newKnownNames=null;
            if(nNames.length==0) break;
        }
    }
    

    /// adds the given input field to the known ones
    void addName(InputField i){
        auto a=i.name in knownNames;
        if (a is null){
            knownNames[i.name]=i;
            newKnownNames~=i.name;
        } else {
            if (!((*a) is i)){
                // duplicated name, allow and add as next name, or check better for equality??
                throw new Exception("duplicated name "~i.name~" in input",__FILE__,__LINE__);
            }
        }
    }
    
    /// returns the unique representant for the field i (the old value, if present, otherwise registers the new one)
    InputField makeUnique(InputField i){
        auto oldV=i.name in knownNames;
        if (oldV !is null){
            oldV.unify(i);
            return *oldV;
        } else {
            knownNames[i.name]=i;
            newKnownNames~=i.name;
            return i;
        }
    }
}

/// reads an xyz file
class FileConfig:Config{
    InputField _myField;
    char[] fileName;
    long frame;
    char[] format;
    ReadSystem readSys;
    ParticleSys!(Real) pSysReal;
    ParticleSys!(LowP) pSysLowP;
    Real[][] cell;
    
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem(){
        if (readSys is null){
            scope dfile=new DataFileInput(fileName);
            scope file=new MultiInput(dfile);
            readSys=readFrame(file,format,frame);
            file.shutdownInput();
            if (cell.length!=0){
                if (cell.length!=3) throw new Exception("invalid cell size (should be 3x3)",__FILE__,__LINE__);
                for (int i=0;i<3;++i) {
                    if (cell[i].length!=3) throw new Exception("invalid cell size (should be 3x3)",__FILE__,__LINE__);
                    for (int j=0;j<3;++j){
                        readSys.cell[i][j]=cell[i][j];
                    }
                }
            }
        }
        return readSys;
    }
    /// the configuration as ParticleSys (might be cached)
    ParticleSys!(Real) particleSysReal(){
        if (pSysReal is null){
            pSysReal=readIn2PSys!(Real)(readSystem());
        }
        return pSysReal;
    }
    ParticleSys!(LowP) particleSysLowP(){
        if (pSysLowP is null){
            pSysLowP=readIn2PSys!(LowP)(readSystem());
        }
        return pSysLowP;
    }
    /// drops the cached readSystem/particleSys
    void clear(){
        readSys=null;
        pSysReal=null;
    }
    mixin myFieldMixin!();
    bool verify(CharSink logger){
        bool res=true;
        auto w=dumper(logger);
        if (fileName.length==0){
            w("Error: fileName is empty in field ")(myFieldName)("\n");
            res=false;
        } else if (!FilePath(fileName).exists()){
            w("Warning: fileName points to a non existing file '")(fileName)("' in field ")(myFieldName)(", continuing expecting it to be created in future\n");
        };
        if (format.length==0){
            if (fileName.length>4){
                switch(fileName[$-4..$]){
                case ".xyz":
                    format="xyz";
                    break;
                case ".car":
                    format="car";
                    break;
                case ".pdb":
                    format="pdb";
                    break;
                default:
                    break;
                }
            }
            if (format.length==0) {
                w("Error: format not given, an could not derive it from fileName in field ")(myFieldName)("\n");
                res=false;
            }
        } else {
            switch(format){
            case "xyz","car","pdb":
                break;
            default:
                w("Error: unknown format '")(format)("' given in field ")(myFieldName)(", known formats are xyz,car,pdb\n");
                res=false;
                break;
            }
        }
        return res;
    }
    mixin(serializeSome("dchem.FileConfig",
    `fileName: the filename where to reat the configuration
    format: the format of the file(xyz,car,pdb)
    frame: if a frame different from the first one should be read (-1 means the last one)
    cell: the cell matrix h (if present overrides the one that might be in the file)`));
    mixin printOut!();
}
