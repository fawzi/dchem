/// calculator
module dchem.calculator.Calculator;
import dchem.Common;
import blip.util.NotificationCenter;
import blip.t.core.Traits: cmp,ctfe_rep;
import blip.t.core.Variant;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;
import blip.parallel.smp.Wait;
import blip.serialization.Serialization;
import blip.sync.UniqueNumber;
import dchem.input.ReadIn;
import dchem.input.RootInput;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import blip.container.BitArray;
import blip.container.Deque;
import dchem.sys.Constraints;
import tango.sys.Process;
import blip.io.IOArray;
import blip.io.NullStream;
import tango.io.vfs.model.Vfs;

/// keeps a list of the allocators for the various methods
class MethodAllocators{
    alias CalculationContext delegate(CalculationInstance cInstance,Method method,char[] className,char[] contextId) MethodAllocator;
    
    MethodAllocator[char[]] methodClasses;
    
    this(){ }
    
    void opIndexAssign(MethodAllocator allocator,char[]name){
        synchronized(this){
            methodClasses[name]=allocator;
        }
    }
    
    MethodAllocator opIndex(char[]name){
        synchronized(this){
            auto r=name in methodClasses;
            if (r is null){
                throw new Exception(collectAppender(delegate void(CharSink s){
                    s("could not find allocator for methodClass '"); s(name); s("'; ");
                    s("known classes:");
                    foreach(k,v;methodClasses){
                        s("'"); s(k); s("'");
                    }
                }),__FILE__,__LINE__);
            }
            return *r;
        }
    }
    
    static MethodAllocators defaultAllocators;
    static this(){
        defaultAllocators=new MethodAllocators();
    }
}

/// represent a manager for instances (execution contextes)
class ClassInstanceManager: Executer {
    bool addToCache;
    CalculationInstance[char[]] activeInstances;
    CalculationInstance[char[]] cache;
    size_t maxInstances=1;
    WaitCondition enoughInstances;
        
    bool verify(CharSink sink){
        bool res=true;
        auto s=dumperP(sink); 
        if (maxInstances==0){
            s("Error: maxInstances should be larget than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    this(size_t maxInstances=1,bool addToCache=true){
        this.maxInstances=maxInstances;
        this.addToCache=addToCache;
        enoughInstances=new WaitCondition(&lessThanMaxInstances);
    }
    bool lessThanMaxInstances(){
        return activeInstances.length<maxInstances;
    }
    CalculationInstance getInstance(InstanceGetFlags flags,uint maxContexts){
        sout("entering ClassInstanceManager get instance\n");
        while (true){
            if ((flags & InstanceGetFlags.ReuseCache)!=0){
                synchronized(this){
                    foreach (instId,inst;cache){
                        auto instance=inst;
                        cache.remove(instId);
                        instance.activate();
                        return instance;
                    }
                }
            }
            if ((flags&InstanceGetFlags.NoAlloc)!=0){
                sout("no alloc\n");
                return null;
            }
            synchronized(this){
                if (activeInstances.length>=maxInstances){
                    if ((flags&InstanceGetFlags.Wait)==0) {
                        sout("no wait\n");
                        return null;
                    }
                } else {
                    sout("create new\n");
                    auto newI=newInstance(maxContexts);
                    newI.activate();
                    return newI;
                }
            }
            sout("begin wait\n");
            enoughInstances.wait();
            sout("try again\n");
        }
    }
    CalculationInstance newInstance(uint maxContexts){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
    }
    
    CalculationInstance getFromCache(char[] instanceId){
        synchronized(this){
            auto i=instanceId in cache;
            if (i !is null){
                auto inst=*i;
                cache.remove(instanceId);
                return inst;
            }
        }
        return null;
    }
    void purgeCache(){
        synchronized(this){
            foreach (k,v;cache){
                v.destroy();
            }
        }
    }
    void activated(CalculationInstance c){
        synchronized(this){
            assert(!(c.instanceId in activeInstances),"activated active entry "~c.instanceId);
            activeInstances[c.instanceId]=c;
        }
    }
    void deactivated(CalculationInstance c){
        bool destroy=false;
        synchronized(this){
            assert(c.instanceId in activeInstances,"deactivated non active entry "~c.instanceId);
            activeInstances.remove(c.instanceId);
            if (addToCache){
                cache[c.instanceId]=c;
            } else {
                destroy=true;
            }
            if (activeInstances.length==maxInstances-1) enoughInstances.notify();
        }
        if (destroy){
            c.destroy();
        }
    }
    void destroyed(CalculationInstance c){
        synchronized(this){
            if (c.instanceId in activeInstances){
                activeInstances.remove(c.instanceId);
                if (activeInstances.length==maxInstances-1) enoughInstances.notify();
            }
            cache.remove(c.instanceId);
        }
    }
    mixin(serializeSome("dchem.Executer",""));
    mixin printOut!();
    mixin myFieldMixin!();
}
/// represent an instance that can calculate systems (i.e. a computational resource)
class CalcInstance:CalculationInstance{
    ClassInstanceManager _manager;
    char[] _instanceId;
    CalculationContext[char[]] contexts; // active contexts
    UniqueNumber!(size_t) lastContextId;
    size_t maxContexts;
    WaitCondition waitForContexts;

    VfsFolder baseDirectory(){
        assert(0,"unimplemented");
        return null;
    }
    
    bool lessThanMaxContexts(){
        return contexts.length<maxContexts;
    }
    this(ClassInstanceManager manager,char[] instanceId,uint maxContexts){
        this._manager=manager;
        this._instanceId=instanceId;
        this.lastContextId=UniqueNumber!(size_t)(1);
        this.maxContexts=maxContexts;
        this.waitForContexts=new WaitCondition(&lessThanMaxContexts);
    }
    Executer manager() { return _manager; }
    char[] instanceId() { return _instanceId; }
    
    CalculationContext newContext(Method m,char[]templateName,CalculationContext delegate(CalculationInstance,Method,char[],char[])allocator){
        char[] contextId=collectAppender(delegate void(CharSink s){
            s(templateName);
            s("-");
            writeOut(s,lastContextId.next());
        });
        while (true){
            synchronized(this){
                if (contexts.length<maxContexts) {
                    auto res=allocator(this,m,templateName,contextId);
                    res.activate();
                    return res;
                }
            }
            waitForContexts.wait();
        }
    }
    void localExec(RemoteTask r){
        r.execute(Variant(this));
    }
    /// creates a process that will be able to execute the given command
    Process cmd(char[] cmd){
        auto p=new Process(cmd);
        p.workDir=baseDirectory.toString();
        return p;
    }
    /// gets output (and errors) of the execution of the given process
    static char[] getOutput(Process p,out int status){
        auto buf=new IOArray(512,512);
        p.redirect=Redirect.Output|Redirect.ErrorToOutput;
        p.execute();
        buf.copy(p.stdout);
        auto res=p.wait();
        status=res.status;
        switch (res.reason)
        {
            case res.Exit:
                break;
            case res.Signal,res.Stop,res.Continue,res.Error:
            default:
                if (status==0) status=-1;
                break;
        }
        p.stdout.close();
        return cast(char[])buf.slice();
    }
    /// executes the given process, discards any output and return the status
    static int execProcess(Process p){
        auto buf=nullStream();
        p.redirect=Redirect.Output|Redirect.ErrorToOutput;
        p.execute();
        buf.copy(p.stdout);
        auto res=p.wait();
        auto status=res.status;
        switch (res.reason)
        {
            case res.Exit:
                break;
            case res.Signal,res.Stop,res.Continue,res.Error:
            default:
                if (status==0) status=-1;
                break;
        }
        p.stdout.close();
        return status;
    }
    void execCmd(char[] cmd,CharSink log=sout.call){
        if (cmd.length>0 && cmd!="NONE"){
            char[256] buf;
            auto arr=lGrowableArray(buf,0);
            dumperP(&arr)("executing ")(cmd)("\n");
            log(arr.data);
            arr.clearData();
            auto p=new Process(cmd);
            int status;
            log(getOutput(p,status));
            if (status!=0){
                arr("command failed with status ");
                writeOut(&arr.appendArr,status);
                throw new Exception(arr.data.dup,__FILE__,__LINE__);
            }
        }
    }
    void activate(){
        manager.activated(this);
    }
    void deactivate(){
        manager.deactivated(this);
    }
    void destroy(){
        manager.destroyed(this);
    }
    void activatedContext(CalculationContext c){
        synchronized(this){
            assert(!(c.contextId in contexts),"activated already active context "~c.contextId);
            contexts[c.contextId]=c;
        }
    }
    void deactivatedContext(CalculationContext c){
        synchronized(this){
            assert((c.contextId in contexts),"deactivated already inactive context "~c.contextId);
            contexts.remove(c.contextId);
            waitForContexts.checkCondition();
        }
    }
    void destroyedContext(CalculationContext c){
        contexts.remove(c.contextId);
        waitForContexts.checkCondition();
    }
    // comparison
    override equals_t opEquals(Object o){
        if (o.classinfo !is this.classinfo){
            return 0;
        }
        auto c=cast(CalcInstance)o;
        if (c !is this){
            assert(c.instanceId!=instanceId,"duplicate instanceId");
            return false;
        }
        return true;
    }
    override hash_t toHash(){
        return getHash(instanceId);
    }
    override int opCmp(Object o){
        if (o.classinfo !is this.classinfo){
            return ((cast(void*)o.classinfo<cast(void*)this.classinfo)?-1:1);
        }
        auto t=cast(CalcInstance)o;
        return cmp(instanceId,t.instanceId);
    }
}

char[] withPSys(char[]op,char[]from=""){
    return `
    if(`~from~`pSysReal !is null){
        auto pSys=`~from~`pSysReal;
        `~op~`
    } else if (`~from~`pSysLowP !is null){
        auto pSys=`~from~`pSysLowP;
        `~op~`
    } else {
        throw new Exception("no valid particle system in context "~`~from~`contextId~" trying to execute "~`~ctfe_rep(op)~`,__FILE__,__LINE__);
    }
    `;
}
/// represent a calculation that might have been aready partially setup, in particular the
/// number of elements,... cannot change
class CalcContext:CalculationContext{
    char[] _contextId;
    ParticleSys!(Real) _pSysReal;
    ParticleSys!(LowP) _pSysLowP;
    SegmentedArray!(Vector!(Real,3)) posArr;
    NotificationCenter _nCenter;
    HistoryManager!(LowP) _posHistory;
    CalculationInstance _cInstance;
    ChangeLevel _changeLevel; /// 0: first time, 1: all changed, 2: only pos changed, 3: small pos change, 4: extrapolable pos change
    Real maxChange;
    MultiConstraint _constraints;
    
    this(CalculationInstance cInstance,char[] contextId){
        this._cInstance=cInstance;
        this._contextId=contextId;
        _nCenter=new NotificationCenter();
        _constraints=new MultiConstraint();
    }
    MultiConstraint constraints(){ return _constraints; }
    char[] contextId(){ return _contextId; }
    ParticleSys!(Real) pSysReal() { return _pSysReal; }
    ParticleSys!(LowP) pSysLowP() { return _pSysLowP; }
    NotificationCenter nCenter() { return _nCenter; }
    HistoryManager!(LowP) posHistory() { return _posHistory; }
    CalculationInstance cInstance() { return _cInstance; }
    ChangeLevel changeLevel() { return _changeLevel; }
    void changeLevel(ChangeLevel c) { _changeLevel=c; }
    Real potentialEnergy(){
        mixin(withPSys("return pSys.dynVars.potentialEnergy;"));
    }
    void pos(SegmentedArray!(Vector!(Real,3)) newPos){
        mixin(withPSys("pSys.dynVars.x.pos[]=newPos;"));
    }
    SegmentedArray!(Vector!(Real,3)) pos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.x.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.x.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.x.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    void dpos(SegmentedArray!(Vector!(Real,3)) newDpos){
        mixin(withPSys("pSys.dynVars.dx.pos[]=newDpos;"));
    }
    SegmentedArray!(Vector!(Real,3)) dpos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.dx.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.dx.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.dx.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    void mddpos(SegmentedArray!(Vector!(Real,3)) newMddpos){
        mixin(withPSys("pSys.dynVars.mddx.pos[]=newMddpos;"));
    }
    SegmentedArray!(Vector!(Real,3)) mddpos(){
        if (pSysReal!is null) {
            return pSysReal.dynVars.mddx.pos;
        } else if (pSysLowP!is null) {
            if (posArr is null){
                posArr=pSysLowP.dynVars.mddx.pos.dupT!(Vector!(Real,3))();
            } else {
                pSysLowP.dynVars.mddx.pos.dupTo(posArr);
            }
            return posArr;
        } else {
            throw new Exception("missing particle sys in context "~contextId,__FILE__,__LINE__);
        }
    }
    SysStruct sysStruct(){
        SysStruct res;
        mixin(withPSys("res=pSys.sysStruct;"));
        return res;
    }
    void changedDynVars(ChangeLevel changeLevel,Real diff){
        if (changeLevel<_changeLevel) _changeLevel=changeLevel;
        if (diff>maxChange) maxChange=diff;
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
        // maxChange=0.0; changeLevel=ChangeLevel.SmoothPosChange;
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

