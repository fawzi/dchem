/// calculator
module dchem.calculator.Calculator;
import dchem.Common;
import blip.util.NotificationCenter;
import blip.core.Traits: cmp,ctfe_rep;
import blip.core.Variant;
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
public import dchem.calculator.CalculatorModels;
import blip.parallel.rpc.Rpc;
import dchem.input.WriteOut;
import blip.parallel.mpi.MpiModels;

/// Limits the number of contexts that can be created/active
class ContextLimiter:InputElement{
    CalculationContext[char[]] contexts; // created contexts
    CalculationContext[char[]] activeContexts; // active contexts
    size_t maxContexts=1;
    size_t maxActiveContexts=size_t.max;
    WaitCondition waitForContexts;
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.ContextLimiter",`maxActiveContexts : maximum number of active contexts (non binding)
    maxContexts : maximum number of live (i.e. not yet given back) contexts`));
    
    bool lessThanMaxContexts(){
        bool res;
        synchronized(this){
            res=contexts.length<maxContexts && activeContexts.length<maxActiveContexts;
        }
        return res;
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
            LocalCalculationContext ctx;
            bool addMore=false;
            synchronized(this){
                if (contexts.length<maxContexts) {
                    addMore=true;
                    --maxContexts;
                }
            }
            if (addMore){
                auto ctx0=cContext(wait,history);
                if (ctx0 is null) return null;
                ctx=cast(LocalCalculationContext)cast(Object)ctx0;
                if (ctx is null) throw new Exception("limiter works only with local contexts",__FILE__,__LINE__);
                synchronized(this){
                    assert((ctx.contextId in contexts)is null,"context already present "~ctx.contextId);
                    assert((ctx.contextId in activeContexts)is null,"context already present "~ctx.contextId);
                    contexts[ctx.contextId]=ctx;
                    activeContexts[ctx.contextId]=ctx;
                    ++maxContexts;
                }
            }
            if(ctx!is null){
                ctx.nCenter.registerCallback("willActivateContext",&willActivateContextCB,Callback.Flags.ReceiveAll);
                ctx.nCenter.registerCallback("willDeactivateContext",&willDeactivateContextCB,Callback.Flags.ReceiveAll);
                ctx.nCenter.registerCallback("willGiveBackContext",&willGiveBackContextCB,Callback.Flags.ReceiveAll);
                return ctx;
            }
            if(!wait){
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
        if (method is null || cast(Method)method.contentObj() is null){
            res=false;
            log("method has to be valid, and a method in ")(myFieldName)("\n");
        }
        return res;
    }
    
    void setup(LinearComm pEnv,CharSink log){ }
    
    /// gets a calculator to perform calculations with this method, if possible reusing the given history
    CalculationContext getCalculator(bool wait, ubyte[]history){
        assert(cl!is null && cast(Method)method.contentObj() !is null);
        auto m=cast(Method)method.contentObj;
        auto dlg=&m.getCalculator;
        return cl.createContext(&m.getCalculator,wait,history);
    }
    /// drops the history with the given id
    void dropHistory(ubyte[]history){
        (cast(Method)method.contentObj).dropHistory(history);
    }
    /// clears all history
    void clearHistory(){
        (cast(Method)method.contentObj).clearHistory();
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
const char[] calcCtxMethodsStr=
`activePrecision|contextId|pSysWriterReal|pSysWriterLowP|pSysWriterRealSet|pSysWriterLowPSet|refPSysReal|refPSysLowP|constraintGen|sysStruct|changeLevel|changeLevelSet|changedDynVars|potentialEnergy|posSet|pos|dposSet|dpos|mddposSet|mddpos|updateEF|activate|deactivate|giveBack|stop|method|storeHistory|exportedUrl|executeLocally`;

/// represent a calculation that might have been aready partially setup, in particular the
/// number of elements,... cannot change
class CalcContext:LocalCalculationContext{
    char[] _contextId;
    ParticleSys!(Real) _pSysReal;
    ParticleSys!(LowP) _pSysLowP;
    SegmentedArray!(Vector!(Real,3)) posArr;
    NotificationCenter _nCenter;
    HistoryManager!(LowP) _posHistory;
    ChangeLevel _changeLevel; /// 0: first time, 1: all changed, 2: only pos changed, 3: small pos change, 4: extrapolable pos change
    Method _method;
    Real maxChange;
    
    void executeLocally(RemoteCCTask r){
        r.workOn(this);
    }

    ConstraintI!(Real) constraintsReal(){ return null; }
    ConstraintI!(LowP) constraintsLowP(){ return null; }
    char[] contextId(){ return _contextId; }
    Precision activePrecision() { return Precision.Real; }
    ParticleSys!(Real) refPSysReal() { return _pSysReal; }
    ParticleSys!(LowP) refPSysLowP() { return _pSysLowP; }
    ParticleSys!(Real) pSysReal() { return _pSysReal; }
    ParticleSys!(LowP) pSysLowP() { return _pSysLowP; }
    PSysWriter!(Real) pSysWriterReal(){ return pSysWriter(_pSysReal); }
    PSysWriter!(LowP) pSysWriterLowP() { return pSysWriter(_pSysLowP); }
    void pSysWriterRealSet(PSysWriter!(Real)p) { assert(_pSysReal!is null); _pSysReal[]=p; }
    void pSysWriterLowPSet(PSysWriter!(LowP)p) { assert(_pSysLowP!is null); _pSysLowP[]=p; }
    ConstraintGen constraintGen(){
        auto cReal=constraintsReal();
        if (cReal!is null){
            return cReal.constraintGen();
        }
        auto cLowP=constraintsLowP();
        if (cLowP !is null){
            return cLowP.constraintGen();
        }
        return null;
    }
    
    NotificationCenter nCenter()  { return _nCenter; }
    HistoryManager!(LowP) posHistory() { return _posHistory; }
    ChangeLevel changeLevel() { return _changeLevel; }
    void changeLevelSet(ChangeLevel c) { _changeLevel=c; }
    Real potentialEnergy(){
        mixin(withPSys("return pSys.dynVars.potentialEnergy;"));
    }
    void posSet(SegmentedArray!(Vector!(Real,3)) newPos){
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
    void dposSet(SegmentedArray!(Vector!(Real,3)) newDpos){
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
    void mddposSet(SegmentedArray!(Vector!(Real,3)) newMddpos){
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
    
    void setup(LinearComm pEnv,CharSink log){ }
    
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
    mixin(rpcMixin("dchem.CalcContext", "CalculationContext",calcCtxMethodsStr));
    DefaultVendor vendor;
    /// url to access this from other processes (only as CalculationContext)
    char[] exportedUrl(){
        return vendor.proxyObjUrl();
    }
    this(char[] contextId){
        this._contextId=contextId;
        _nCenter=new NotificationCenter();
        // register to the world...
        vendor=new DefaultVendor(this);
        assert(ProtocolHandler.defaultProtocol!is null,"defaultProtocol");
        assert(ProtocolHandler.defaultProtocol.publisher!is null,"publisher");
        ProtocolHandler.defaultProtocol.publisher.publishObject(vendor,"CalcContext"~contextId,true,Publisher.Flags.Public);
    }
}

