/// calculator
module dchem.calculator.Calculator;
import blip.util.NotificationCenter;

interface RemoteTask:Serializable{
    void execute(Variant args);
}

enum InstanceGetFlags{
    ReuseCache=1,    /// reuse instances that are in the cache
    NoAlloc=2,       /// don't alloc a new instance
    NoAllocSubOpt=4, /// don't alloc a new instance if suboptimal
    Wait=5           /// wait until an instance is available (when one hits the maximum number of instances)
}

class CalcInstanceManager{
    ClassInstanceManager[char[]] instanceManagers;
    
    void registerClass(char[]name,ClassInstanceManager cManager){
        synchronize(this){
            instanceManagers[name]=cManager;
        }
    }
    
    CalcInstance getInstanceForClass(char[]iClass,InstanceGetFlags flags){
        synchronize(this){
            return instanceManagers[name].getInstanceForClass;
        }
    }
}

/// represent a manager for instances
class ClassInstanceManager{
    this(){}
    CalcInstance getInstanceForClass(char[]iClass,InstanceGetFlags flags){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
    }
    CalcInstance getFromCache(char[] instanceId){
        return null;
    }
    void putInCache(CalcInstance cInstance){
        if (cInstance !is null){
            cInstance.destroy();
        }
    }
    void purgeCache(){}
}
/// represent an instance that can calculate systems
class CalcInstance{
    char[] instanceId;
    char[] instanceClass;

    this(char[] instanceId,char[]instanceClass){
        this.instanceId=instanceId;
        this.instanceClass=instanceClass;
    }
    CalcContext newContext(sysStruct,char[]templateName){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
    }
    CalcContext getFromCache(char[] contextId){
        return null;
    }
    void putInCache(CalcContext context){
        context.destroy;
    }
    void purgeCache(){}
    bool hasLocalExecution(){ return false; }
    void localExec(RemoteTask r){
        throw new Exception("remote execution not supported",__FILE__,__LINE__);
    }
    void destroy(){}
}

/// keeps a history of the previous steps
class HistoryManager{
    size_t nHistory; /// number of steps to keep
    BitArray keepH; /// things to keep in history
    enum bPos{ /// bit positions in keepH
        cell=0,
        pos=1,
        dpos=2,
        ddpos=3
        orient=4,
        dorient=5,
        ddorient=6,
        dof=7,
        ddof=8,
        dddof=9
    }
    Deque!(ParticleSys) history; /// place to keep the history
    this(size_t nHistory=1,BitArray keepH=BitArray([true,true,false,false,true,false,false,true,false,false])){
        this.nHistory=nHistory;
        this.keepH=keepH;
        history=new Deque!(ParticleSys)(nHistory);
    }
    void addToHistory(ParticleSys p){
        if (history.length<nHistory){
            DynamicsVars dVars;
            if (keepH[bPos.cell])
                dVars.cell=p.cell.dup;
            if (keepH[bPos.pos])
                dVars.pos=p.pos.dup;
            if (keepH[bPos.dpos])
                dVars.dpos=p.dpos.dup;
            if (keepH[bPos.ddpos])
                dVars.ddpos=p.ddpos.dup;
            if (keepH[bPos.orient])
                dVars.orient=p.orient.dup;
            if (keepH[bPos.dorient])
                dVars.dorient=p.dorient.dup;
            if (keepH[bPos.ddorient])
                dVars.ddorient=p.ddorient.dup;
            if (keepH[bPos.dof])
                dVars.dof=p.dof.dup;
            if (keepH[bPos.ddof])
                dVars.ddof=p.ddof.dup;
            if (keepH[bPos.dddof])
                dVars.dddof=p.dddof.dup;
            ParticleSys n=new ParticleSys(collectAppender(delegate(CharSink s){
                s("history-"); writeOut(history.length); }),
                p.sysStruct,dVars,cast(NotificationCenter)null);
            history.pushFront(n);
        } else {
            auto nO=history.popEnd();
            nO[]=p;
            history.pushFront(nO);
        }
    }
}
/// represent a calculation that has been aready partially setup, in particular the
/// number of elements,... cannot change
class CalcContext{
    char[] contextId;
    ParticleSys pSys;
    NotificationCentral nCentral;
    CalcInstance cInstance;
    Real maxChange;
    int changeLevel; /// 0: first time, 1: all changed, 2: only pos changed, 3: small pos change, 4: extrapolable pos change
    
    this(char[] contextId,SysStruct sysStruct,CalcInstance cInstance){
        this.contextId=contextId;
        this.sysStruct=sysStruct;
        this.cInstance=cInstance;
        nCentral=new NotificationCentral();
    }
    
    void pos(SegmentedArray!(Real) newPos){
        pSys.dynVars.pos[]=newPos;
    }
    SegmentedArray!(Real) pos(){
        return pSys.dynVars.pos;
    }
    void dpos(SegmentedArray!(Real) newDpos){
        pSys.dynVars.dpos[]=newDpos;
    }
    SegmentedArray!(Real) dpos(){
        return pSys.dynVars.dpos;
    }
    void pos(SegmentedArray!(Real) newDdpos){
        pSys.dynVars.ddpos[]=newDdpos;
    }
    SegmentedArray!(Real) ddpos(){
        return pSys.dynVars.ddpos;
    }
    
    void dynVars(DynamicsVars *newDv){
        return pSys.dynVars[]=*newDv;
    }
    DynamicsVars *dynVars(){
        return pSys.dynVars;
    }
    
    void ChangedDynVars(int changeLevel,Real diff){
        
    }
    
/+  // use dynVars instead of these...
    // orientation (quaternions)
    SegmentedArray!(Real[4]) orient(){
        return pSys.dynVars.orient;
    }
    SegmentedArray!(Real[4]) dorient(){
        pSys.dynVars.dorient;
    }
    SegmentedArray!(Real[4]) ddorient(){
        pSys.dynVars.ddorient;
    }
    void orient(SegmentedArray!(Real[4])newVal){
        pSys.dynVars.orient[]=newVal;
    }
    void dorient(SegmentedArray!(Real[4])newVal){
        pSys.dynVars.dorient[]=newVal;
    }
    void ddorient(SegmentedArray!(Real[4])newVal){
        pSys.dynVars.ddorient[]=newVal;
    }

    // other degrees of freedom to be integrated (internal coords, extended lagrangian,...)
    SegmentedArray!(Real) dof(){
        return pSys.dynVars.dof;
    }
    SegmentedArray!(Real) ddof(){
        return pSys.dynVars.ddof;
    }
    SegmentedArray!(Real) dddof(){
        return pSys.dynVars.dddof;
    }
    void dof  (SegmentedArray!(Real) newVal){
        pSys.dynVars.dof[]=newVal;
    }
    void ddof (SegmentedArray!(Real) newVal){
        pSys.dynVars.ddof[]=newVal;
    }
    void dddof(SegmentedArray!(Real) newVal){
        pSys.dynVars.dddof[]=newVal;
    }
    +/
    
    Real e_pot(){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
    }
    
    void updateEF(bool updateE=true,bool updateF=true){
        throw new Exception("to implement in subclasses",__FILE__,__LINE__);
        // maxChange=0.0; changeLevel=0;
    }
    
    void destroy(){}
}