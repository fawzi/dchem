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
import dchem.sys.DynVars;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import blip.container.Deque;
import blip.container.BulkArray;
import dchem.sys.Constraints;
import dchem.sys.PIndexes;
import tango.sys.Process;
import blip.io.IOArray;
import blip.io.NullStream;
import tango.io.vfs.model.Vfs;
public import dchem.calculator.CalculatorModels;
import blip.parallel.rpc.Rpc;
import dchem.input.WriteOut;
import blip.parallel.mpi.MpiModels;
import PosUtils=dchem.sys.PosUtils;
import dchem.sys.Cell;
import blip.sync.Atomic;

/// cell based distances for positions
///
/// an object that handles the various distance measures, and distance related things
/// it defines the macroscopic distance space, the microscopic distance space (for small changes)
/// is defined in ParticleSys by its overlap, here one defines things like PBC.
class NDistOps:DistOpsGen,DistOps{
    bool useFirstNeigh=true;
    mixin(serializeSome("dchem.NDistOps",`
    useFirstNeigh: if the first image convention should be used for periodic directions (true)`));
    mixin printOut!();
    mixin myFieldMixin!();
    this(char[] fName=""){
        if (fName.length>0){
            myField=new InputField(fName,InputField.TypeId.InputElement,this);
        }
    }
    bool verify(CharSink s){
        return true;
    }
    
    DistOps DistOpsForContext(CalculationContext){
        return this;
    }
    
    /// wraps deltaX (a difference between two X points) so that it is as small as possible
    /// compatibly with the implicit symmetries (but not with the explicit ones)
    void wrapT(T)(ParticleSys!(T)pSys,DynPVector!(T,XType)deltaX){
        PosUtils.wrap(pSys.dynVars.x.cell,deltaX.pos);
    }
    // alias don't work reliably :(
    // ditto
    void wrap(ParticleSys!(Real)pSys,DynPVector!(Real,XType)deltaX){
        wrapT!(Real)(pSys,deltaX);
    }
    /// ditto
    void wrap(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)deltaX){
        wrapT!(LowP)(pSys,deltaX);
    }
    /// ditto
    void wrapReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)deltaX){
        wrapT!(Real)(pSys,deltaX);
    }
    /// ditto
    void wrapLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)deltaX){
        wrapT!(LowP)(pSys,deltaX);
    }

    /// get the periodic copy of x that is closest (in first image convetion, i.e. in the h_inv distance)
    /// to pSys (this might be different than the closest one for very skewed h)
    void makeCloseT(T)(ParticleSys!(T)pSys,DynPVector!(T,XType)x){
        PosUtils.makeClose(pSys.dynVars.x.cell,pSys.dynVars.x.pos,x.pos);
    }
    void makeClose(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x){
        makeCloseT!(Real)(pSys,x);
    }
    /// ditto
    void makeClose(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x){
        makeCloseT!(LowP)(pSys,x);
    }
    /// ditto
    void makeCloseReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x){
        makeCloseT!(Real)(pSys,x);
    }
    /// ditto
    void makeCloseLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x){
        makeCloseT!(LowP)(pSys,x);
    }

    /// distance in the reduced units, this is the norm2 distance between x and pSys.dynVars.x after a makeClose
    /// or deltaX.norm2 after calling wrap, where deltaX=x-pSys.dynVars.x
    Real reducedDistT(T)(ParticleSys!(T)pSys,DynPVector!(T,XType)x,Real threshold){
        Real totDist=0;
        if (useFirstNeigh){
            auto dAtt=PosUtils.configDist(pSys.dynVars.x.cell,pSys.dynVars.x.pos,x.pos);
            atomicAdd(totDist,dAtt);
        } else {
            auto o1=pSys.dynVars.x.pos;
            auto o2=x.pos;
            assert(o1.kRange==o2.kRange);
            foreach(k;o1.kRange.pLoop){
                auto dAtt=PosUtils.norm22Threshold(o1[k].basicData(),o2[k].basicData(),threshold-totDist);
                atomicAdd(totDist,dAtt);
            }
        }
        if (totDist<threshold){
            auto o1=pSys.dynVars.x.orient;
            auto o2=x.orient;
            assert(o1.kRange==o2.kRange);
            foreach(k;o1.kRange.pLoop){
                auto dAtt=PosUtils.norm22Threshold(o1[k].basicData(),o2[k].basicData(),threshold-totDist);
                atomicAdd(totDist,dAtt);
            }
        }
        if (totDist<threshold){
            auto d1=pSys.dynVars.x.dof;
            auto d2=x.dof;
            assert(d1.kRange==d2.kRange);
            foreach(k;d1.kRange.pLoop){
                auto dAtt=PosUtils.norm22Threshold(d1[k].basicData(),d2[k].basicData(),threshold-totDist);
                atomicAdd(totDist,dAtt);
            }
        }
        return totDist;
    }
    Real reducedDist(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold){
        return reducedDistT!(Real)(pSys,x,threshold);
    }
    /// ditto
    Real reducedDist(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold){
        return reducedDistT!(LowP)(pSys,x,threshold);
    }
    /// ditto
    Real reducedDistReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold){
        return reducedDistT!(Real)(pSys,x,threshold);
    }
    /// ditto
    Real reducedDistLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold){
        return reducedDistT!(LowP)(pSys,x,threshold);
    }

    /// distance in reduced units between *one* particle of kinf k at pos1,ord1,dof1, and possibly many
    /// particles of the same kind at pos2,ord2,dof2
    void rDistOneToNT(T)(ParticleSys!(T)pSys,KindIdx k,
        BulkArray!(Vector!(T,3))pos1,BulkArray!(Quaternion!(T))ord1,BulkArray!(T)dof1,
        BulkArray!(Vector!(T,3))pos2,BulkArray!(Quaternion!(T))ord2,BulkArray!(T)dof2,
        T[] dists)
    {
        auto cell=pSys.dynVars.x.cell;
        PosUtils.oneToNDistKind!(T)(cell.h.cell,cell.hInv.cell,cell.periodicFlags,pos1.basicData(),
            pos2.basicData(),dists,cell.flags);
        PosUtils.addNorm22OneToN(ord1.basicData,ord2.basicData,dists);
        PosUtils.addNorm22OneToN(dof1.basicData,dof2.basicData,dists);
    }
    void rDistOneToN(ParticleSys!(Real)pSys,KindIdx k,
        BulkArray!(Vector!(Real,3))pos1,BulkArray!(Quaternion!(Real))ord1,BulkArray!(Real)dof1,
        BulkArray!(Vector!(Real,3))pos2,BulkArray!(Quaternion!(Real))ord2,BulkArray!(Real)dof2,
        Real[] dists){
        rDistOneToNT!(Real)(pSys,k,pos1,ord1,dof1,pos2,ord2,dof2,dists);
    }
    /// ditto
    void rDistOneToN(ParticleSys!(LowP)pSys,KindIdx k,
        BulkArray!(Vector!(LowP,3))pos1,BulkArray!(Quaternion!(LowP))ord1,BulkArray!(LowP)dof1,
        BulkArray!(Vector!(LowP,3))pos2,BulkArray!(Quaternion!(LowP))ord2,BulkArray!(LowP)dof2,
        LowP[] dists){
        rDistOneToNT!(LowP)(pSys,k,pos1,ord1,dof1,pos2,ord2,dof2,dists);
    }
    /// ditto
    void rDistOneToNReal(ParticleSys!(Real)pSys,KindIdx k,
        BulkArray!(Vector!(Real,3))pos1,BulkArray!(Quaternion!(Real))ord1,BulkArray!(Real)dof1,
        BulkArray!(Vector!(Real,3))pos2,BulkArray!(Quaternion!(Real))ord2,BulkArray!(Real)dof2,
        Real[] dists){
        rDistOneToNT!(Real)(pSys,k,pos1,ord1,dof1,pos2,ord2,dof2,dists);
    }
    /// ditto
    void rDistOneToNLowP(ParticleSys!(LowP)pSys,KindIdx k,
        BulkArray!(Vector!(LowP,3))pos1,BulkArray!(Quaternion!(LowP))ord1,BulkArray!(LowP)dof1,
        BulkArray!(Vector!(LowP,3))pos2,BulkArray!(Quaternion!(LowP))ord2,BulkArray!(LowP)dof2,
        LowP[] dists){
        rDistOneToNT!(LowP)(pSys,k,pos1,ord1,dof1,pos2,ord2,dof2,dists);
    }

    /// full (cartesian) distance, this is the reducedDist in the full system using cartesian units
    /// in general it might be expensive to calculate, or be not better than the 
    /// if threshold is different from 0, then as soon as the distance is detected to be larger than it
    /// it is returned, which might be quicker if one does not care when for the exact value of distances
    /// larger than threshold
    Real fullDist(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x,Real threshold=0){
        return reducedDist(pSys,x,threshold);
    }
    /// ditto
    Real fullDist(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x,Real threshold=0){
        return reducedDist(pSys,x,threshold);
    }
    /// ditto
    Real fullDistReal(ParticleSys!(Real)pSys,DynPVector!(Real,XType)x, Real threshold=0){
        return reducedDist(pSys,x,threshold);
    }
    /// ditto
    Real fullDistLowP(ParticleSys!(LowP)pSys,DynPVector!(LowP,XType)x, Real threshold=0){
        return reducedDist(pSys,x,threshold);
    }
}

/// an object that loops on all symmetry equivalent structures generated by neigh that are within epsilon of pSys
/// when having no symmetry
class NoSymmNeighLooper:SymmNeighLooperGen,SymmNeighLooper{
    mixin(serializeSome("dchem.NoSymmNeighLooper",``));
    mixin myFieldMixin!();
    mixin printOut!();
    this(string fName=""){
        if (fName.length>0) {
            myField=new InputField(fName,InputField.TypeId.InputElement,this);
        }
    }
    bool verify(CharSink s){
        return true;
    }
    SymmNeighLooper symmNeighLooperForContext(CalculationContext c){
        return this;
    }
    int loopOnNeighWithinT(T)(ParticleSys!(T)pSys,DistOps distOps,DynPVector!(T,XType)neigh,T epsilon,
        int delegate(ref DynPVector!(T,XType))loopBody)
    {
        auto dist=distOps.reducedDist(pSys,neigh,epsilon);
        if (dist<epsilon){
            return loopBody(neigh);
        }
    }
    // aliases don't work reliably
    int loopOnNeighWithin(ParticleSys!(Real)pSys,DistOps distOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody)
    {
        return loopOnNeighWithinT!(Real)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithin(ParticleSys!(LowP)pSys,DistOps distOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody)
    {
        return loopOnNeighWithinT!(LowP)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithinReal(ParticleSys!(Real)pSys,DistOps distOps,DynPVector!(Real,XType)neigh,Real epsilon,
        int delegate(ref DynPVector!(Real,XType))loopBody)
    {
        return loopOnNeighWithinT!(Real)(pSys,distOps,neigh,epsilon,loopBody);
    }
    int loopOnNeighWithinLowP(ParticleSys!(LowP)pSys,DistOps distOps,DynPVector!(LowP,XType)neigh,LowP epsilon,
        int delegate(ref DynPVector!(LowP,XType))loopBody)
    {
        return loopOnNeighWithinT!(LowP)(pSys,distOps,neigh,epsilon,loopBody);
    }
}

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
    mixin printOut!();
    
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
    CharSink _log;
    CharSink logger(){
        return _log;
    }
    
    this(){
        _log=sout.call;
    }
    
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.ContextLimiterClient",`contextLimiter : the limiter that constrain this method
    method : the method that is constrained when accessed through this`));
    mixin printOut!();

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
    
    void setup(LinearComm pEnv,CharSink log){
        _log=log;
        if (method !is null){
            method.contentT!(Method)().setup(pEnv,log);
        }
    }
    
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
    switch(`~from~`activePrecision){
    case Precision.Real:
        auto pSys=`~from~`pSysReal;
        `~op~`
        break;
    case Precision.LowP:
        auto pSys=`~from~`pSysLowP;
        `~op~`
        break;
    default:
        throw new Exception("unexpected precision "~`~from~`contextId~" trying to execute "~`~ctfe_rep(op)~`,__FILE__,__LINE__);
    }
    `;
}
const char[] calcCtxMethodsStr=
`activePrecision|contextId|pSysWriterReal|pSysWriterLowP|pSysWriterRealSet|pSysWriterLowPSet|`~
`refPSysReal|refPSysLowP|constraintGen|sysStruct|changeLevel|changeLevelSet|changedDynVars|`~
`potentialEnergy|potentialEnergyError|posSet|pos|dposSet|dpos|mddposSet|mddpos|updateEF|activate|deactivate|`~
`giveBack|stop|method|storeHistory|exportedUrl|executeLocally|logMsg1`;

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
    CharSink _logger;
    SymmNeighLooper _symmNeighLooper;
    DistOps _distOps;
    
    void logMsg1(char[]m){
        _logger(m);
    }
    CharSink logger(){
        return _logger;
    }
    
    void executeLocally(RemoteCCTask r){
        r.workOn(this);
    }
    
    void desc(CharSink s){
        dumper(s)("executer '")(contextId())("'");
    }

    ConstraintI!(Real) constraintsReal(){ return null; }
    ConstraintI!(LowP) constraintsLowP(){ return null; }
    SymmNeighLooper symmNeighLooper(){
        if (_symmNeighLooper is null){
            synchronized(this){
                if (_symmNeighLooper is null){
                    _symmNeighLooper=new NoSymmNeighLooper("defaultNoSymmLooper");
                }
            }
        }
        return _symmNeighLooper;
    }
    DistOps distOps(){
        if (_distOps is null){
            synchronized(this){
                if (_distOps is null){
                    _distOps=new NDistOps("defaultNDistOps");
                }
            }
        }
        return _distOps;
    }
    char[] contextId(){ return _contextId; }
    Precision activePrecision() { return Precision.Real; }
    
    ParticleSys!(Real) refPSysReal() {
        if (_pSysReal is null){
            assert(activePrecision==Precision.LowP);
            auto newPSys=_pSysLowP.dupT!(Real)(PSDupLevel.EmptyDyn);
            synchronized(this){
                if (_pSysReal is null){
                    _pSysReal=newPSys;
                } else {
                    newPSys.release();
                }
            }
        }
        return _pSysReal;
    }
    ParticleSys!(LowP) refPSysLowP() {
        if (_pSysLowP is null){
            assert(activePrecision==Precision.Real);
            auto newPSys=_pSysReal.dupT!(LowP)(PSDupLevel.EmptyDyn);
            synchronized(this){
                if (_pSysLowP is null){
                    _pSysLowP=newPSys;
                } else {
                    newPSys.release();
                }
            }
        }
        return _pSysLowP;
    }
    void sysStructChangedCallback(cstring n,Callback*cb,Variant v){
        assert(n=="sysStructChanged","unexpected callback name:"~n);
        switch(activePrecision){
        case Precision.Real:
            if (_pSysLowP is null){
                _pSysLowP.release();
            }
            _pSysLowP=null;
            break;
        case Precision.LowP:
            if (_pSysReal is null){
                _pSysReal.release();
            }
            _pSysReal=null;
            break;
        default:
            assert(0);
        }
    }
    ParticleSys!(Real) pSysReal() { if (activePrecision==Precision.Real) return _pSysReal; return null; }
    ParticleSys!(LowP) pSysLowP() { if (activePrecision==Precision.LowP) return _pSysLowP; return null; }
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
    Real potentialEnergyError(){
        mixin(withPSys("return pSys.dynVars.potentialEnergyError;"));
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
    this(char[] contextId,CharSink log,NotificationCenter nCenter=null){
        this._contextId=contextId;
        this._logger=log;
        this._nCenter=nCenter;
        if (nCenter is null)
            this._nCenter=new NotificationCenter();
        this._nCenter.registerCallback("sysStructChanged",&this.sysStructChangedCallback,Callback.Flags.ReceiveAll);
        // register to the world...
        vendor=new DefaultVendor(this);
        assert(ProtocolHandler.defaultProtocol!is null,"defaultProtocol");
        assert(ProtocolHandler.defaultProtocol.publisher!is null,"publisher");
        ProtocolHandler.defaultProtocol.publisher.publishObject(vendor,"CalcContext"~contextId,true,Publisher.Flags.Public);
    }
}

