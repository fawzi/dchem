/// Describes an object that can loop on neighbors (neighboring list)
module dchem.neigh.NeighModels;
import blip.BasicModels;
import blip.serialization.Serialization;
import blip.math.IEEE;
import dchem.Common;
import dchem.sys.Cell;
import dchem.sys.PIndexes;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;

alias float SafeScalar;
alias Vector!(float,3) SafeVector;
alias Quaternion!(float) SafeQuaternion;

template SafeType(T){
    static if (is(T==float)||is(T==double)||is(T==real)){
        alias float SafeType;
    } else static if (is(T==Vector!(T.flt,T.dim))){
        alias Vector!(SafeType!(T.flt),T.dim) SafeType;
    } else static if (is(T==Quaternion!(T.flt))){
        alias Quaternion!(SafeType!(T.flt)) SafeType;
    } else static if (is(T==NArray!(T.dtype,T.dim))){
        alias NArray!(SafeType!(T.dtype),T.dim) SafeType;
    } else {
        static assert(0,"type "~T~" not supported by SafeType");
    }
}

SafeType!(T) toSafeType(T)(T a){
    static if (is(T==SafeType!(T))){
        return a;
    } else static if (is(T==float)||is(T==double)||is(T==real)){
        return cast(float)a;
    } else static if (is(T==Vector!(T.flt,T.dim))){
        return SafeType!(T).from(a);
    } else static if (is(T==Quaternion!(T.flt))){
        return SafeType!(T).from(a);
    } else static if (is(T==NArray!(T.dtype,T.dim))){
        return a.asType!(SafeType!(T.dtype))();
    } else {
        static assert(0,"type "~T~" not supported by toSafeType");
    }
}

/// defines cutoffs for neighboring lists
/// this is the optimized struct that is used for loops, only uses KindIdx
/// if you build it programmatically you probably want to use this
struct DCutoffs{
    Real _maxR;
    Real[KindIdx] kindIdxR;
    Real defaultR;
    bool upperTriangle=true;
    
    Real maxR(){
        if (isNaN(_maxR)){
            _maxR=defaultR;
            foreach(k,v;kindIdxR){
                if (!(v<=_maxR)){
                    _maxR=v;
                }
            }
        }
        return _maxR;
    }
    Real rForKindIdx(KindIdx k){
        auto r1=k in kindIdxR;
        if (r1!is null){
            return *r1;
        }
        return defaultR;
    }
    mixin(serializeSome("dchem.DCutoffs","Cutoff radii for the different kinds","kindIdxR|defaultR"));
    mixin printOut!();
}

/// defines cutoffs for neighboring lists
/// this is the "main" struct that can be defined for example in the input,
/// as it accepts for example string labels.
struct DistCutoffs{
    Real _maxR;
    Real[string] kindR;
    Real[string] symbolR;
    Real[KindIdx] kindIdxR;
    Real defaultR;
    bool upperTriangle=true;
    /// restricts the cutoffs to the given SysStruct (discarding unneeded )
    DCutoffs restrictTo(SysStruct sys){
        DCutoffs res;
        foreach(k,v;symbolR){
            foreach(kind;sys.kinds){ // avoid quadratic scaling?
                if (kind.symbol==k){
                    res.kindIdxR[kind.pKind]=v;
                    break;
                }
            }
        }
        foreach(k,v;kindR){
            foreach(kind;sys.kinds){ // avoid quadratic scaling?
                if (kind.name==k){
                    res.kindIdxR[kind.pKind]=v;
                    break;
                }
            }
        }
        foreach(k,v;kindIdxR){
            res.kindIdxR[k]=v;
        }
        res.defaultR=defaultR;
        res.upperTriangle=upperTriangle;
        return res;
    }
    
    Real maxR(){
        if (isNaN(_maxR)){
            _maxR=defaultR;
            foreach(k,v;kindR){
                if (!(v<=_maxR)){
                    _maxR=v;
                }
            }
            foreach(k,v;symbolR){
                if (!(v<=_maxR)){
                    _maxR=v;
                }
            }
            foreach(k,v;kindIdxR){
                if (!(v<=_maxR)){
                    _maxR=v;
                }
            }
        }
        return _maxR;
    }
    
    Real rKind(ParticleKind pK){
        auto r1=pK.name in kindR;
        if (r1 !is null){
            return *r1;
        }
        auto r2=pK.symbol in symbolR;
        if (r2!is null){
            return *r2;
        }
        return defaultR;
    }
    mixin(serializeSome("dchem.DistCutoffs","Cutoff radii for the different kinds","kindR|symbolR|kindIdxR|defaultR"));
    mixin printOut!();
}

/// generates neighboring lists
interface NeighListGen:BasicObjectI{
    /// creates the structure to build neighboring lists on the given configuration 
    NeighList!(Real) createOnConfigReal(Cell!(Real) cell,SegmentedArray!(Vector!(Real,3)) sArr);
    /// ditto
    NeighList!(LowP) createOnConfigLowP(Cell!(LowP) cell,SegmentedArray!(Vector!(LowP,3)) sArr);
    /// ditto
    NeighList!(Real) createOnConfig(Cell!(Real) cell,SegmentedArray!(Vector!(Real,3)) sArr);
    /// ditto
    NeighList!(LowP) createOnConfig(Cell!(LowP) cell,SegmentedArray!(Vector!(LowP,3)) sArr);
}

/// describes a pair of neighbors
struct NeighPair(T){
    PIndex p1;
    index_type i1;
    PIndex p2;
    index_type i2;
    Vector!(T,3) p1p2;
}

/// a neighboring list
/// parallel loops *might* be parallel, but don't have to
interface NeighList(T):BasicObjectI{
    alias SafeType!(T) SafeT;
    enum :size_t { defaultOptimalBlockSize=((16*1024<T.sizeof)?cast(size_t)1:16*1024/T.sizeof) }
    
    /// range of kinds this neighboring list loops on
    KindRange kRange();
    
    alias int delegate(ref NeighPair!(T) neighs) LoopBody;
    
    /// loops on the neighbor of p with kind k
    int sloopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2,LoopBody loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize);
    /// loops on the neighbors between particles of kind k1 and k2
    int sloopOnNeighsKK(KindIdx k1,KindIdx k2,SafeT cutoff2,LoopBody loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize);
    
    /// parallel loop on the neighbor of p with kind k
    int ploopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2,LoopBody loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize);
    /// parallel loop on the neighbors between particles of kind k1 and k2
    int ploopOnNeighsKK(KindIdx k1,KindIdx k2,SafeT cutoff2, LoopBody loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize);
    
    /// reatin/relase semantic
    void retain();
    /// ditto
    void release();
}

/// helper structure to loop on neighbors
struct NeighLooper(T){
    NeighList!(T) looper;
    enum :size_t { defaultOptimalBlockSize=NeighList!(T).defaultOptimalBlockSize }
    
    static struct LoopPK(loopType){
        NeighList!(T) looper;
        SafeT cutoff2;
        PIndex p1;
        index_type i1;
        KindIdx k;
        bool avoidDoubleCount=false;
        size_t optSize=defaultOptimalBlockSize;
        
        int opApply(int delegate(ref NeighPair!(T) neighs) loopBody){
            static if ((loopType&1)!=0) {
                looper.ploopOnNeighsPK(p1,i1,k,cutoff2,loopBody,avoidDoubleCount,optSize);
            } else {
                looper.sloopOnNeighsPK(p1,i1,k,cutoff2,loopBody,avoidDoubleCount,optSize);
            }
        }
        static if ((loopType&1)==0) {
            int opApply(int delegate(ref ulong, ref NeighPair!(T)) loopBody){
                ulong idx=0;
                looper.sloopOnNeighsPK(p1,i1,k,cutoff2,delegate int(ref NeighPair!(T)neighs){
                    int res=loopBody(idx,neighs);
                    if (res) return res;
                    ++idx;
                },avoidDoubleCount,optSize);
            }
        }
    }
    
    /// loops on the neighbors of p with kind k
    LoopPK!(loopType) loopOnNeighsPK(int loopType)(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        LoopPK!(loopType) res;
        res.cutoff2=cutoff2;
        res.p1=p1;
        res.i1=i1;
        res.k=k;
        res.avoidDoubleCount=avoidDoubleCount;
        res.optSize=optSize;
        res.looper=looper;
        return res;
    }
    /// sequential loop on the neighbors of p with kind k
    LoopPK!(LoopType.Sequential) sloopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2, bool avoidDoubleCount=false)
    {
        return loopOnNeighsPK!(LoopType.Sequential)(p1,i1,k,cutoff2,avoidDoubleCount);
    }
    /// parallel loop on the neighbors of p with kind k
    LoopPK!(LoopType.Parallel) ploopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return loopOnNeighsPK!(LoopType.Parallel)(p1,i1,k,cutoff2,avoidDoubleCount,optSize);
    }
    
    
    static struct LoopP(loopType){
        NeighList!(T) looper;
        PIndex p1;
        index_type i1;
        bool avoidDoubleCount=false;
        size_t optSize=defaultOptimalBlockSize;
        DCutoffs *dCutoffs;
        
        int opApply(int delegate(ref NeighPair!(T) neighs) loopBody){
            SafeT cut1=dCutoffs.rForKindIdx(p1.kind);
            static if ((loopType&1)!=0) {
                foreach (k; looper.kRange) {
                    looper.ploopOnNeighsPK(p1,i1,k,pow2(cut1+dCutoffs.rForKindIdx(k)),loopBody,avoidDoubleCount,optSize);
                }
            } else {
                foreach (k; looper.kRange.pLoop) {
                    looper.sloopOnNeighsP(p1,i1,k,pow2(cut1+dCutoffs.rForKindIdx(k)),loopBody,avoidDoubleCount);
                }
            }
        }
        static if ((loopType&1)==0) {
            int opApply(int delegate(ref ulong, ref NeighPair!(T)) loopBody){
                ulong idx=0;
                SafeT cut1=dCutoffs.rForKindIdx(p1.kind);
                foreach (k; looper.kRange) {
                    looper.sloopOnNeighsP(p1,i1,k,pow2(cut1+dCutoffs.rForKindIdx(k)),delegate int(ref NeighPair!(T)neighs){
                        int res=loopBody(idx,neighs);
                        if (res) return res;
                        ++idx;
                    },avoidDoubleCount,optSize);
                }
            }
        }
    }
    
    /// loops on the neighbors of p
    LoopP!(loopType) loopOnNeighsP(int loopType)(PIndex p1,index_type i1, DCutoffs *dCutoffs, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        LoopP!(loopType) res;
        res.p1=p1;
        res.i1=i1;
        res.avoidDoubleCount=avoidDoubleCount;
        res.optSize=optSize;
        res.looper=looper;
        res.dCutoffs=dCutoffs;
        return res;
    }
    /// sequential loop on the neighbors of p
    LoopP!(LoopType.Sequential) sloopOnNeighsP(PIndex p1,index_type i1,DCutoffs *dCutoffs, bool avoidDoubleCount=false)
    {
        return loopOnNeighsP!(LoopType.Sequential)(p1,i1,dCutoffs,avoidDoubleCount);
    }
    /// parallel loop on the neighbors of p
    LoopP!(LoopType.Parallel) ploopOnNeighsP(PIndex p1,index_type i1, DCutoffs *dCutoffs ,bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return loopOnNeighsP!(LoopType.Parallel)(p1,i1,dCutoffs,avoidDoubleCount,optSize);
    }
    
    static struct LoopKK(loopType){
        NeighList!(T) looper;
        KindIdx k1;
        KindIdx k2;
        SafeT cutoff2;
        bool avoidDoubleCount=false;
        size_t optSize=defaultOptimalBlockSize;
        
        int opApply(int delegate(ref NeighPair!(T) neighs) loopBody){
            static if ((loopType&1)!=0) {
                looper.ploopOnNeighsKK(k1,k2,cutoff2,loopBody,avoidDoubleCount,optSize);
            } else {
                looper.sloopOnNeighsKK(k1,k2,cutoff2,loopBody,avoidDoubleCount,optSize);
            }
        }
        static if ((loopType&1)==0) {
            int opApply(int delegate(ref ulong, ref NeighPair!(T)) loopBody){
                ulong idx=0;
                looper.sloopOnNeighsKK(k1,k2,cutoff2,delegate int(ref NeighPair!(T)neighs){
                    int res=loopBody(idx,neighs);
                    if (res) return res;
                    ++idx;
                },avoidDoubleCount,optSize);
            }
        }
    }
    
    /// loops on the neighbor of particles of kind k1 that have kind k2
    LoopKK!(loopType) loopOnNeighsKK(int loopType)(KindIdx k1, KindIdx k2, SafeT cutoff2, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        LoopKK!(loopType) res;
        res.k1=k1;
        res.k2=k2;
        res.cutoff2=cutoff2;
        res.avoidDoubleCount=avoidDoubleCount;
        res.optSize=optSize;
        res.looper=looper;
        return res;
    }
    /// sequential loop on the neighbors of p with kind k
    LoopKK!(LoopType.Sequential) sloopOnNeighsKK(KindIdx k1, KindIdx k2, bool avoidDoubleCount=false)
    {
        return loopOnNeighsKK!(LoopType.Sequential)(k1,k2,avoidDoubleCount);
    }
    /// parallel loop on the neighbors of p with kind k
    LoopKK!(LoopType.Parallel) ploopOnNeighsKK(KindIdx k1, KindIdx k2, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return loopOnNeighsKK!(LoopType.Parallel)(k1,k2,avoidDoubleCount,optSize);
    }
    
    static struct LoopNeighs(loopType){
        NeighList!(T) looper;
        bool avoidDoubleCount=false;
        size_t optSize=defaultOptimalBlockSize;
        DCutoffs *dCutoffs;
        
        int opApply(int delegate(ref NeighPair!(T) neighs) loopBody){
            static if ((loopType&1)!=0) {
                foreach(k1;loopType.kRange.ploop){
                    SafeT r1=dCutoffs.rForKindIdx(k1);
                    foreach(k2;loopType.kRange.ploop){
                        SafeT r2=dCutoffs.rForKindIdx(k2);
                        looper.ploopOnNeighsKK(k1,k2,pow2(r1+r2),loopBody,avoidDoubleCount,optSize);
                    }
                }
            } else {
                foreach(k1;loopType.kRange){
                    SafeT r1=dCutoffs.rForKindIdx(k1);
                    foreach(k2;loopType.kRange){
                        SafeT r2=dCutoffs.rForKindIdx(k2);
                        looper.sloopOnNeighs(loopBody,avoidDoubleCount,optSize);
                    }
                }
            }
        }
        static if ((loopType&1)==0) {
            int opApply(int delegate(ref ulong, ref NeighPair!(T)) loopBody){
                ulong idx=0;
                foreach(k1;loopType.kRange){
                    SafeT r1=dCutoffs.rForKindIdx(k1);
                    foreach(k2;loopType.kRange){
                        SafeT r2=dCutoffs.rForKindIdx(k2);
                        looper.sloopOnNeighs(k1,k2,pow2(r1+r2),delegate int(ref NeighPair!(T)neighs){
                                int res=loopBody(idx,neighs);
                                if (res) return res;
                                ++idx;
                            },avoidDoubleCount,optSize);
                    }
                }
            }
        }
    }
    
    /// loops on all the neighbors
    LoopNeighs!(loopType) loopOnNeighs(int loopType)(bool avoidDoubleCount=false, size_t optSize=defaultOptimalBlockSize)
    {
        LoopNeighs!(loopType) res;
        res.avoidDoubleCount=avoidDoubleCount;
        res.optSize=optSize;
        res.looper=looper;
        return res;
    }
    /// sequential loop on the neighbors of p with kind k
    LoopNeighs!(LoopType.Sequential) sloopOnNeighs(bool avoidDoubleCount=false)
    {
        return loopOnNeighs!(LoopType.Sequential)(avoidDoubleCount);
    }
    /// parallel loop on the neighbors of p with kind k
    Loop!(LoopType.Parallel) ploopOnNeighs(bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return loopOnNeighs!(LoopType.Parallel)(avoidDoubleCount,optSize);
    }
}