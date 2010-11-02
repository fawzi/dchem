/// describes the root input file
module dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.text.TextParser;
import blip.io.BasicIO;
import blip.io.Console;
import dchem.Common;
import blip.util.TangoLog;
import blip.core.Variant;
import blip.BasicModels;
import blip.parallel.mpi.MpiModels;

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

/// a calculation method
interface Sampler:InputElement{
    void run(LinearComm pWorld,CharSink log);
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
        Sampler,
        Reference,
        InputField,
        InputElement,
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
    
    Object contentObj(){
        return cast(Object)content;
    }
    T contentT(T)(bool raiseOnNull=true){
        auto obj=cast(T)cast(Object)content;
        if (obj is null && raiseOnNull){
            throw new Exception("content of "~myFieldName~" cannot be casted to "~T.stringof,__FILE__,__LINE__);
        }
        return obj;
    }

    static char[]typeStr(TypeId t){
        switch (t){
            case TypeId.Sampler: return "sampler";
            case TypeId.Reference: return "ref";
            case TypeId.InputField: return "inputField";
            case TypeId.InputElement: return "inputElement";
            default: assert(0,"invalid typeId");
        }
    }
    static TypeId strType(char[] val){
        switch (val){
            case "sampler": return TypeId.Sampler;
            break;
            case "ref":     return TypeId.Reference;
            break;
            case "inputField": return TypeId.InputField;
            break;
            case "inputElement": return TypeId.InputElement;
            break;
            default: assert(0,"invalid typeId");
        }
    }
    /// sets typeId depending on content
    static TypeId objType (Serializable o){
        if (o is null) return TypeId.Reference;
        auto oo=cast(Object)o;
        if (cast(Sampler)oo) return TypeId.Sampler;
        if (cast(InputField)oo) return TypeId.InputField;
        if (cast(InputElement)oo) return TypeId.InputElement;
        throw new Exception("could not find type of object",__FILE__,__LINE__);
    }
    Sampler sampler(){
        return cast(Sampler)cast(Object)content;
    }
    InputField inputField(){
        return cast(InputField)cast(Object)content;
    }
    InputElement inputElement(){
        return cast(InputElement)cast(Object)content;
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


/// generates a helper tmplated function to generate type baseType~I!(T) from an (input generated)
/// baseType~Gen generator (a common pattern for the input objects that have pecision dependent
/// template specializations)
char[] genTypeTMixin(char[]baseType,char[]object,char[] extraArgsDecl,char[] extraArgs){
    char[] res=`
    `~baseType~`I!(T) `~object~`T(T)(`~baseType~`Gen gen`~((extraArgsDecl.length>0)?",":" ")~extraArgsDecl~`){
        static if (is(T==Real)){
            return gen.`~object~`Real(`~extraArgs~`);
        } else static if (is(T==LowP)){
            return gen.`~object~`LowP(`~extraArgs~`);
        } else {
            static assert(0,T.stringof~" is neither Real nor LowP");
        }
    }`;
    return res;
}
