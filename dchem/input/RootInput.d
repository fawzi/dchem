/// describes the root input file
module dchem.input.RootInput;
import dchem.sys.CoreDefs:ParticleSys;

/// configuration
interface Config:Serializable{
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem();
    /// the configuration as ParticleSys (might be cached)
    ParticleSys particleSys();
    /// drops the cached readSystem/particleSys
    void clear();
}

/// calculator setup
interface Calculator:Serializable{
    void setupCalculatorClass();
    char[] calculatorClass();
}

/// a calculation method
interface Method:Serializable{
    void run();
    void stop();
}

/// an root input entry
/// serializes as reference, unserializes normally
class InputField:Serializable{
    enum TypeId{
        Configuration,
        Method,
        Sampler,
        Reference
    }
    char[] name;
    Serializable content;
    TypeId typeId;
    InputField next;

    char[]typeStr(){
        switch (typeId){
            case TypeId.Configuration: return "config";
            case TypeId.Method: return "method";
            case TypeId.Sampler: return "sampler";
            case TypeId.Reference: return "ref";
            default: assert(0,"invalid typeId");
        }
    }
    void typeStr(char[] val){
        switch (typeId){
            case "config":  typeId=TypeId.Configuration;
            break;
            case "method":  typeId=TypeId.Method;
            break;
            case "sampler": typeId=TypeId.Sampler;
            break;
            case "ref":     typeId=TypeId.Reference;
            break;
            default: assert(0,"invalid typeId");
        }
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
            Log.lookup("blip.parallel.smp").warning("warning: serializing InputField without a RootInput in serializer context",
                __FILE__,__LINE__);
            rInp=new RootInput();
            s.context["RootInput"]=rInp;
        } else {
            rInp=rInp0.get!(RootInput)();
        }
        rInp.addName(name,this);
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
            Log.lookup("blip.parallel.smp").warning("warning: unserializing InputField without a RootInput in serializer context",
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
                content=o.name;
            } else {
                if (o.typeId!=typeId){
                    throw new Exception("unification of InputFields of different types",__FILE__,__LINE__);
                }
                if (content!is o.content && !overwrite){
                    throw new Exception("unification of InputFields changed value and overwrite is false",__FILE__,__LINE__);
                }
                content=o.content;
            }
        }
    }
    /// really serializes the content of this 
    struct RealSerial{
        InputField ptr;
        static ClassMetaInfo metaI;
        static this(){
            metaI=ClassMetaInfo.createForType!(typeof(*this))("InputField");// use T.mangleof?
            metaI.kind=TypeKind.CustomK;
        }
        ClassMetaInfo getSerializationMetaInfo(){
            return metaI;
        }
        static if (is(typeof(ptr.preSerialize(Serializer.init)))){
            void preSerialize(Serializer s){
                if (ptr) ptr.preSerialize(s);
            }
        }
        static if (is(typeof(ptr.postSerialize(Serializer.init)))){
            void postSerialize(Serializer s){
                if (ptr) ptr.postSerialize(s);
            }
        }
        void serialize(Serializer s){
            FieldMetaInfo *elMetaInfoP=null;
            version(PseudoFieldMetaInfo){
                FieldMetaInfo elMetaInfo=FieldMetaInfo("el","",
                    getSerializationInfoForType!(InputField)());
                elMetaInfo.pseudo=true;
                elMetaInfoP=&elMetaInfo;
            }
            if (ptr is null){
                s.field(elMetaInfoP,null);
            } else {
                s.field(elMetaInfoP, ptr.content);
            }
        }
        static if (is(typeof(ptr.preUnserialize(Unserializer.init)))){
            Serializable preUnserialize(Unserializer s){
                if (ptr) return ptr.preUnserialize(s);
                return this;
            }
        }
        /// unserialize an object
        void unserialize(Unserializer s){
            FieldMetaInfo *elMetaInfoP=null;
            version(PseudoFieldMetaInfo){
                FieldMetaInfo elMetaInfo=FieldMetaInfo("el","",
                    getSerializationInfoForType!(InputField)());
                elMetaInfo.pseudo=true;
                elMetaInfoP=&elMetaInfo;
            }
            if (ptr is null){
                s.field(elMetaInfoP,null);
            } else {
                s.field(elMetaInfoP, ptr.content);
            }
        }
        static if (is(typeof(ptr.postUnserialize(Unserializer.init)))){
            Serializable postUnserialize(Unserializer s){
                if (ptr) return ptr.postUnserialize(s);
                return this;
            }
        }
    }
}


class RootInput{
    InputField[char[]] knownNames;
    InputField[char[]] newKnownNames;
    bool rewriteNamesOk;
    
    /// reads the input from the given TextParser
    void readInput(TextParser!(char) t){
        scope js=new JsonUnserializer!(char)(t)
        js.context["RootInput"]=Variant(this);
        char[] name;
        js.handlers(name);
        if (t.getSeparator()!=":"){
            throw new Exception("unexpected separator",__FILE__,__LINE__);
        }
        auto inputF=new InputField(name);
        js(inputF.content);
        rootNames[name]=inputF;
        if (inputF.content!is null){
            solvedRef(name,inputF.content);
        }
    }

    /// adds the given input field to the known ones
    void addName(InputField i){
        auto a=i.name in knownNames;
        if (a is null){
            knownNames[i.name]=i;
            newKnownNames[i.name]=i;
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
            return oldV;
        } else {
            knownNames[i.name]=i;
            return i;
        }
    }
}

/// reads an xyz file
class FileConfig:Config{
    char[] fileName;
    long frame;
    char[] format;
    ReadSystem readSys;
    ParticleSys pSys;
    
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem(){
        if (readSys is null){
            switch(format){
            case "xyz":
                auto file=new TextParser!(char)(new DataFileInput(fileName));
                file.newlineIsSpace=false;
                for (long iframe=1;iframe<frame;++iframe){
                    size_t nat;
                    file(nat);
                    tp.skipLines(nat+2);
                }
                readSys=readXYZFrame(file);
                if (frame==-1){
                    while (true){
                        auto readSys1=readXYZFrame(file,readSys,true);
                        if (readSys1 is null) break;
                        readSys=readSys1;
                    }
                }
                break;
            case "car":
                // fully disallow empty frames???
                readSys=readCarHeader(file);
                auto readSys1=readCarFrame(file,readSys,true);
                if (readSys1 !is null) readSys=readSys1;
                for (long iframe=1;iframe<frame;++iframe){
                    readSys1=readCarFrame(file,readSys,true);
                    if (readSys1 !is null) readSys=readSys1;
                }
                if (frame==-1){
                    while (true){
                        readSys1=readXYZFrame(file,readSys,true);
                        if (readSys1 is null) break;
                        readSys=readSys1;
                    }
                }
                break;
            case "pdb":
                if (frame!=0) throw new Exception("only single frame pdb are supported",__FILE__,__LINE__);
                readSys=readPdb(file,readSys,true);
                break;
            default:
                assert(0,"unexpected fromat "~format);
            }
        }
        return readSys;
    }
    /// the configuration as ParticleSys (might be cached)
    ParticleSys particleSys();
    /// drops the cached readSystem/particleSys
    void clear();
}