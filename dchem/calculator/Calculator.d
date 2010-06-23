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
import blip.container.Deque;
import dchem.sys.Constraints;
import tango.sys.Process;
import blip.io.IOArray;
import blip.io.NullStream;
import tango.io.vfs.model.Vfs;

/// Limits the number of contexts that can be created/active
class ContextLimiter:InputElement{
    CalculationContext[char[]] contexts; // created contexts
    CalculationContext[char[]] activeContexts; // active contexts
    size_t maxContexts;
    size_t maxActiveContexts;
    WaitCondition waitForContexts;
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.ContextLimiter",`maxActiveContexts : maximum number of active contexts (non binding)
    maxContexts : maximum number of live (i.e. not yet given back) contexts`));
    
    bool lessThanMaxContexts(){
        return contexts.length<maxContexts && activeContexts.length<maxActiveContexts;
    }
    bool verify(CharSink sink){
        auto s=dumper(sink);
        bool res=true;
        if (maxActiveContexts<=0){
            s("maxActiveContexts must be larger than 0 in ")(myFieldName)("\n");
            res=false;
        }
        if (maxContexts<=0){
            s("maxContexts must be larger than 0 in ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    this(size_t maxContexts=1,size_t maxActiveContexts=1){
        assert(maxContexts!=0,"maxContexts cannot be 0");
        assert(maxActiveContexts!=0,"maxActiveContexts cannot be 0");
        this.maxContexts=maxContexts;
        this.maxActiveContexts=maxActiveContexts;
        this.waitForContexts=new WaitCondition(&lessThanMaxContexts);
    }
    
    /// creates a context if there aren't too many, if wait is true waits until some can be created
    CalculationContext createContext(CalculationContext delegate(bool,ubyte[]history)cContext,
        bool wait,ubyte[]history){
        while (true){
            CalculationContext ctx;
            synchronized(this){
                if (contexts.length<maxContexts) {
                    ctx=cContext(wait,history);
                    if(ctx!is null){
                        contexts[ctx.contextId]=ctx;
                        activeContexts[ctx.contextId]=ctx;
                        ctx.nCenter.registerCallback("willActivateContext",&willActivateContextCB,Callback.Flags.ReceiveAll);
                        ctx.nCenter.registerCallback("willDeactivateContext",&willDeactivateContextCB,Callback.Flags.ReceiveAll);
                        ctx.nCenter.registerCallback("willGiveBackContext",&willGiveBackContextCB,Callback.Flags.ReceiveAll);
                    }
                    return ctx;
                }
            }
            if(ctx!is null || (!wait)){
                return ctx;
            }
            waitForContexts.wait();
        }
    }

    void willActivateContextCB(char[]name,Callback*callBack,Variant v){
        auto c=v.get!(CalculationContext)();
        synchronized(this){
            assert((c.contextId in contexts),"activated context not in list "~c.contextId);
            assert((c.contextId in activeContexts),"activated context already activated "~c.contextId);
            activeContexts[c.contextId]=c;
        }
    }
    void willDeactivateContextCB(char[]name,Callback*callBack,Variant v){
        auto c=v.get!(CalculationContext)();
        synchronized(this){
            assert((c.contextId in activeContexts),"deactivated context not in list "~c.contextId);
            activeContexts.remove(c.contextId);
        }
        waitForContexts.checkCondition();
    }
    void willGiveBackContextCB(char[]name,Callback*callBack,Variant v){
        auto c=v.get!(CalculationContext)();
        synchronized(this){
            assert((c.contextId in contexts),"gave back context not in list "~c.contextId);
            if((c.contextId in activeContexts)!is null)
                activeContexts.remove(c.contextId);
            contexts.remove(c.contextId);
        }
        waitForContexts.checkCondition();
    }
}

/// wraps a method so that it is constrained to a maximum number of contexts when used though this method
class ContextLimiterClient:Method{
    InputField contextLimiter;
    InputField method;
    
    this(){
        // should expose it to the world...
    }
    
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.ContextLimiterClient",`contextLimiter : the limiter that constrain this method
    method : the method that is constrained when accessed through this`));

    final ContextLimiter cl(){
        return cast(ContextLimiter)cast(Object)contextLimiter.content;
    }

    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        if (contextLimiter is null || cl() is null){
            res=false;
            log("contextLimiter has to be valid, and of type dchem.ContextLimiter in ")(myFieldName)("\n");
        }
        if (method is null || method.method() is null){
            res=false;
            log("method has to be valid, and a method in ")(myFieldName)("\n");
        }
        return res;
    }
    /// gets a calculator to perform calculations with this method, if possible reusing the given history
    CalculationContext getCalculator(bool wait, ubyte[]history){
        assert(cl!is null && method.method!is null);
        return cl.createContext(&method.method.getCalculator,wait,history);
    }
    /// drops the history with the given id
    void dropHistory(ubyte[]history){
        method.method.dropHistory(history);
    }
    /// clears all history
    void clearHistory(){
        method.method.clearHistory();
    }
    /// url to access this from other processes
    char[] exportedUrl(){
        assert(0,"to do");
        return "";
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
    ChangeLevel _changeLevel; /// 0: first time, 1: all changed, 2: only pos changed, 3: small pos change, 4: extrapolable pos change
    Method _method;
    Real maxChange;
    MultiConstraint _constraints;
    
    void localExec(RemoteTask r){
        r.execute(Variant(this));
    }

    this(char[] contextId){
        this._contextId=contextId;
        _nCenter=new NotificationCenter();
        _constraints=new MultiConstraint();
        // register to the world...
    }
    MultiConstraint constraints(){ return _constraints; }
    char[] contextId(){ return _contextId; }
    ParticleSys!(Real) pSysReal() { return _pSysReal; }
    ParticleSys!(LowP) pSysLowP() { return _pSysLowP; }
    NotificationCenter nCenter()  { return _nCenter; }
    HistoryManager!(LowP) posHistory() { return _posHistory; }
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
    
    /// called automatically after creation, but before any energy evaluation
    /// should be called before working again with a deactivated calculator
    void activate(){
        nCenter.notify("willActivateContext",Variant(this));
    }
    /// call this to possibly get rid of all caches (i.e. before a pause in the calculation)
    void deactivate(){
        nCenter.notify("willDeactivateContext",Variant(this));
    }
    /// call this to remove the context (after all calculations with this are finished)
    void giveBack(){
        nCenter.notify("willGiveBackContext",Variant(this));
    }
    /// tries to stop a calculation in progress. Recovery after this is not possible
    /// giveBack should still be called
    void stop(){}
    /// the method of this calculator (this might be different from the method that was used to create this
    /// as that might have been wrapped)
    Method method(){
        assert(0,"to implement in subclasses");
    }
    /// stores the history somewhere and returns an id to possibly recover that history at a later point
    /// this is just an optimization, it does not have to do anything. If implemented then
    /// method.getCalculator, .dropHistory and .clearHistory have to be implemented accordingly
    ubyte[]storeHistory(){ return []; }
    /// exposes (publish/vends) this object to the world
    void publish(){
    }
    /// url to access this from other processes
    char[] exportedUrl(){
        assert(0,"to do");
    }
}

