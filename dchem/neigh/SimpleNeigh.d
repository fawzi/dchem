/// loops on the  neighbors using the simple O(N**2) algorithm
module dchem.neigh.SimpleNeigh;
import dchem.Common;
import dchem.neigh.NeighModels;
import dchem.sys.Cell;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import dchem.sys.PIndexes;
import blip.io.BasicIO;
import blip.io.Console;
import blip.util.RefCount;
import blip.math.IEEE;
import blip.math.Math;
import blip.container.BulkArray;
import blip.container.GrowableArray;
import blip.serialization.Serialization;
import blip.util.Convert;
import blip.narray.NArray;
//r12 = (r1_l-r2_l+h_lm*I_m)*(r1_l-r2_l+h_lm*I_m) =
//  = (r1_l-r2_l)**2+2*(r1_l-r2_l)*h_lm*I_m+h_ln*h_lk*I_n*I_k
int loopEllipsoidBoundsT(T, int periodicFlags, bool doTests=false)(Cell!(T) cell, Vector!(T,3) r1, Vector!(T,3) r2, T cutoff,
    PIndex p1, index_type i1, PIndex p2, index_type i2, int delegate(ref NeighPair!(T)) loop){
    alias SafeType!(T) SafeT;
    alias Matrix!(SafeT,3,3) M;
    static if (periodicFlags==0){
        if (pow2(cast(SafeT)(r2.x-r1.x))+pow2(cast(SafeT)(r2.y-r1.y))+pow2(cast(SafeT)(r2.z-r1.z))< cutoff) {
            NeighPair!(T) neigh;
            neigh.p1=p1;
            neigh.i1=i1;
            neigh.p2=p2;
            neigh.i2=i2;
            neigh.p1p2=r2-r1;
            return loop(neigh);
        }
        return 0;
    } else {
        M h = convertTo!(M)(cell.h);
        SafeT mii = h[0,0]*h[0,0]+h[1,0]*h[1,0]+h[2,0]*h[2,0],
            mij = h[0,0]*h[0,1]+h[1,0]*h[1,1]+h[2,0]*h[2,1],
            mik = h[0,0]*h[0,2]+h[1,0]*h[1,2]+h[2,0]*h[2,2],
            mjj = h[0,1]*h[0,1]+h[1,1]*h[1,1]+h[2,1]*h[2,1],
            mjk = h[0,1]*h[0,2]+h[1,1]*h[1,2]+h[2,1]*h[2,2],
            mkk = h[0,2]*h[0,2]+h[1,2]*h[1,2]+h[2,2]*h[2,2];
        SafeT r12_x = cast(SafeT)(r2.x-r1.x);
        SafeT r12_y = cast(SafeT)(r2.y-r1.y);
        SafeT r12_z = cast(SafeT)(r2.z-r1.z);
        SafeT li = r12_x*h[0,0]+r12_y*h[1,0]+r12_z*h[2,0];
        SafeT lj = r12_x*h[0,1]+r12_y*h[1,1]+r12_z*h[2,1];
        SafeT lk = r12_x*h[0,2]+r12_y*h[1,2]+r12_z*h[2,2];
        SafeT a = r12_x*r12_x+r12_y*r12_y+r12_z*r12_z-cutoff;
        // p:=2*li*i+2*lj*j+2*lk*k+i*(mii*i+mij*j+mik*k)+j*(mij*i+mjj*j+mjk*k)+k*(mik*i+mjk*j+mkk*k)+a
        // dp/di = 2*li+2*mii*i+2*mij*j+2*mik*k
        // dp/dj = 2*lj+2*mij*i+2*mjj*j+2*mjk*k
        // dp/dk = 2*lk+2*mik*i+2*mjk*j+2*mkk*k
        assert(periodicFlags==cell.periodicFlags,"unxpected value of periodicFlags");
        int iMin,iMax;
        static if ((periodicFlags & CellPeriodic.x)==0){
            iMin=0;
            iMax=0;
        } else {{
            // eq3 <~> p=0, eq2 <~> dp/dk=0
            SafeT eq3ii=mii, eq3i=2*li, eq3=a;
            SafeT eq3ik=2*mik, eq3k=2*lk, eq3kk=mkk;
            SafeT eq2=lk, eq2i=mik, eq2k=mkk;
            static if ((periodicFlags & CellPeriodic.y)!=0) {
                // substitute dp/dj=0 <=> j=(-lj-mij*i-mjk*k)/mjj
                SafeT fact2 = mjk/mjj;
                eq2  -= fact2*lj;
                eq2i -= fact2*mij;
                eq2k -= fact2*mjk;
                SafeT ljn = lj/mjj;
                eq3  -= ljn*lj;
                eq3i -= 2*ljn*mij;
                eq3k -= 2*ljn*mjk;
                eq3ik -= 2*mij*mjk/mjj;
                eq3kk -= mjk*mjk/mjj;
                eq3ii -= mij*mij/mjj;
            }
            static if ((periodicFlags & CellPeriodic.z)!=0) {
                SafeT valK = -eq2/eq2k;
                SafeT valKI = -eq2i/eq2k;
                eq3   += (eq3kk*valK+eq3k)*valK;
                eq3i  += (2*eq3kk*valK+eq3k)*valKI+eq3ik*valK;
                eq3ii += (eq3kk*valKI+eq3ik)*valKI;
            }
            SafeT delta = eq3i*eq3i-4*eq3ii*eq3;
            if (eq3ii==0 || delta<=0) return 0;
            SafeT sDelta = sqrt(delta);
            static if (doTests) {{
                // test
                SafeT iVal1 = (-eq3i-sDelta)/(2*eq3ii);
                SafeT iVal2 = (-eq3i+sDelta)/(2*eq3ii);
                SafeT kVal1 = 0;
                SafeT kVal2 = 0;
                if ((periodicFlags & CellPeriodic.z)!=0) {
                    kVal1 = -(eq2+eq2i*iVal1)/eq2k;
                    kVal2 = -(eq2+eq2i*iVal2)/eq2k;
                }
                SafeT jVal1 = 0;
                SafeT jVal2 = 0;
                if ((periodicFlags & CellPeriodic.y)!=0) {
                    jVal1 = -(lj+mij*iVal1+mjk*kVal1)/mjj;
                    jVal2 = -(lj+mij*iVal2+mjk*kVal2)/mjj;
                }
                auto iVect1=Vector!(T,3)(iVal1,jVal1,kVal1);
                auto iVect2=Vector!(T,3)(iVal2,jVal2,kVal2);
                Vector!(T,3) r12_1 = r2-r1+cell.h*iVect1;
                Vector!(T,3) r12_2 = r2-r1+cell.h*iVect2;
                if (feqrel2(r12_1.norm22(),cutoff)<SafeT.mant_dig*3/4){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("r12_1 i too far: ")(r12_1.norm22())(" vs ")(cutoff)(" diff:")(r12_1.norm22()-cutoff)
                            (" cell:")(cell)(" r1:")(r1)(" r2:")(r2)(" i:")(iVal1)(" ")(" j:")(jVal1)(" ")(" k:")(kVal1);
                    }));
                }
                if (feqrel2(r12_2.norm22(),cutoff)<SafeT.mant_dig*3/4){
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("r12_2 i too far: ")(r12_2.norm22())(" vs ")(cutoff)(" diff:")(r12_2.norm22()-cutoff)
                            (" cell:")(cell)(" r1")(r1)(" r2:")(r2)(" i:")(iVal2)(" ")(" j:")(jVal2)(" ")(" k:")(kVal2);
                    }));
                }
            }}
            iMin = cast(int)ceil((-eq3i-sDelta)/(2*eq3ii));
            iMax = cast(int)floor((-eq3i+sDelta)/(2*eq3ii));
        }}
        if (iMin>iMax) return 0;
        NeighPair!(T) neigh;
        neigh.p1=p1;
        neigh.i1=i1;
        neigh.p2=p2;
        neigh.i2=i2;
        for (int i=iMin; i<=iMax; ++i) {
            Vector!(T,3) r12i=r2-r1-cell.h.column(0)*iMin;
            int jMin,jMax;
            static if ((periodicFlags & CellPeriodic.y)==0) {
                jMin=0;
                jMax=0;
            } else {{
                SafeT eq3   = a+(2*li+mii*i)*i;
                SafeT eq3j  = 2*lj+2*i*mij;
                SafeT eq3jj = mjj;
                static if ((periodicFlags & CellPeriodic.z)!=0){
                    SafeT valK = -(lk+mik*i)/mkk;
                    SafeT valKj = -mjk/mkk;
                    eq3   += (2*(lk+mik*i)+mkk*valK)*valK;
                    eq3j  += 2*(mik*i+lk+mkk*valK)*valKj+2*mjk*valK;
                    eq3jj += mkk*valKj*valKj+2*mjk*valKj;
                }
                SafeT delta = eq3j*eq3j-4*eq3jj*eq3;
                if (eq3jj==0 || delta<=0) continue;
                SafeT sDelta = sqrt(delta);
                static if (doTests) {{
                    // test
                    SafeT jVal1 = (-eq3j-sDelta)/(2*eq3jj);
                    SafeT jVal2 = (-eq3j+sDelta)/(2*eq3jj);
                    SafeT kVal1 = 0;
                    SafeT kVal2 = 0;
                    static if ((periodicFlags & CellPeriodic.z)!=0) {
                        kVal1 = -(lk+mik*i)/mkk-jVal1*mjk/mkk;
                        kVal2 = -(lk+mik*i)/mkk-jVal2*mjk/mkk;
                    }
                    auto iVect1=Vector!(T,3)(i,jVal1,kVal1);
                    auto iVect2=Vector!(T,3)(i,jVal2,kVal2);
                    Vector!(T,3) r12_1 = r2-r1+cell.h*iVect1;
                    Vector!(T,3) r12_2 = r2-r1+cell.h*iVect2;
                    if (feqrel2(r12_1.norm22(),cutoff)<SafeT.mant_dig*3/4){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("r12_1 j too far: ")(r12_1.norm22())(" vs ")(cutoff)(" diff:")(r12_1.norm22()-cutoff)
                                (" cell:")(cell)(" r1:")(r1)(" r2:")(r2)(" i:")(i)(" ")(" j:")(jVal1)(" ")(" k:")(kVal1);
                        }));
                    }
                    if (feqrel2(r12_2.norm22(),cutoff)<SafeT.mant_dig*3/4){
                        throw new Exception(collectAppender(delegate void(CharSink s){
                            dumper(s)("r12_2 j too far: ")(r12_2.norm22())(" vs ")(cutoff)(" diff:")(r12_2.norm22()-cutoff)
                                (" cell:")(cell)(" r1:")(r1)(" r2:")(r2)(" i:")(i)(" ")(" j:")(jVal2)(" ")(" k:")(kVal2);
                        }));
                    }
                }}
                jMin = cast(int)ceil((-eq3j-sDelta)/(2*eq3jj));
                jMax = cast(int)floor((-eq3j+sDelta)/(2*eq3jj));
            }}
            for (int j=jMin; j<=jMax; ++j) {
                Vector!(T,3) r12j=r12i+cell.h.column(1)*jMin;
                int kMin=0, kMax=0;
                static if ((periodicFlags & CellPeriodic.z)!=0) {{
                    SafeT eq3   = a+2*(li*i+lj*j+mij*i*j)+mii*i*i+mjj*j*j;
                    SafeT eq3k  = 2*(lk+mik*i+mjk*j);
                    SafeT eq3kk = mkk;
                    SafeT delta = eq3k*eq3k-4*eq3kk*eq3;
                    if (eq3kk==0 || delta<=0) continue;
                    SafeT sDelta = sqrt(delta);
                    static if (doTests) {
                        // test 
                        SafeT kVal1 = (-eq3k-sDelta)/(2*eq3kk);
                        SafeT kVal2 = (-eq3k+sDelta)/(2*eq3kk);
                        auto iVect1=Vector!(T,3)(i,j,kVal1);
                        auto iVect2=Vector!(T,3)(i,j,kVal2);
                        Vector!(T,3) r12_1 = r2-r1+h*iVect1;
                        Vector!(T,3) r12_2 = r2-r1+h*iVect2;
                        if (feqrel2(r12_1.norm22(),cutoff)<SafeT.mant_dig*3/4){
                            throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("r12_1 k too far: ")(r12_1.norm22())(" vs ")(cutoff)(" diff:")(r12_1.norm22()-cutoff)
                                    (" cell:")(cell)(" r1:")(r1)(" r2:")(r2)(" i:")(i)(" ")(" j:")(j)(" ")(" k:")(kVal1);
                            }));
                        }
                        if (feqrel2(r12_2.norm22(),cutoff)<SafeT.mant_dig*3/4){
                            throw new Exception(collectAppender(delegate void(CharSink s){
                                dumper(s)("r12_2 k too far: ")(r12_2.norm22())(" vs ")(cutoff)(" diff:")(r12_2.norm22()-cutoff)
                                    (" cell:")(cell)(" r1:")(r1)(" r2:")(r2)(" i:")(i)(" ")(" j:")(j)(" ")(" k:")(kVal2);
                            }));
                        }
                    }
                    kMin = cast(int)ceil((-eq3k-sDelta)/(2*eq3kk));
                    kMax = cast(int)floor((-eq3k+sDelta)/(2*eq3kk));
                }}
                if (kMin>kMax) continue;
                neigh.p1p2=r12j+kMin*cell.h.column(2);
                for (int k=kMin; true; ++k) {
                    int lVal=loop(neigh);
                    if (lVal != 0) return lVal;
                    if (k>kMax) break;
                    neigh.p1p2+=cell.h.column(2);
                }
            }
        }
        return 0;
    }
}

int loopEllipsoidBounds(T, bool doTests=false)(Cell!(T) cell, Vector!(T,3) r1, Vector!(T,3) r2, T cutoff,
    PIndex p1, index_type i1, PIndex p2, index_type i2, int delegate(ref NeighPair!(T)) loop)
{
    switch(cell.periodicFlags){
    case 0:
        return loopEllipsoidBoundsT!(T, 0, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 1:
        return loopEllipsoidBoundsT!(T, 1, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 2:
        return loopEllipsoidBoundsT!(T, 2, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 3:
        return loopEllipsoidBoundsT!(T, 3, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 4:
        return loopEllipsoidBoundsT!(T, 4, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 5:
        return loopEllipsoidBoundsT!(T, 5, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 6:
        return loopEllipsoidBoundsT!(T, 6, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    case 7:
        return loopEllipsoidBoundsT!(T, 7, doTests)(cell, r1, r2, cutoff, p1, i1, p2, i2, loop);
    default:
        assert(0,"unexpected value for periodicFlags in loopEllipsoidBounds");
    }
}
/// loops on the  neighbors using the simple O(N**2) algorithm
class SimpleNeighListGen:NeighListGen{
    this(){}
    /// setup operations to build neighboring lists on the given configuration 
    NeighList!(Real) createOnConfigReal(Cell!(Real) cell, SegmentedArray!(Vector!(Real,3)) sArray){
        return new SimpleNeighList!(Real)(cell, sArray);
    }
    /// ditto
    NeighList!(LowP) createOnConfigLowP(Cell!(LowP) cell, SegmentedArray!(Vector!(LowP,3)) sArray){
        return new SimpleNeighList!(LowP)(cell, sArray);
    }
    /// ditto 
    NeighList!(Real) createOnConfig(Cell!(Real) cell, SegmentedArray!(Vector!(Real,3)) sArray){
        return new SimpleNeighList!(Real)(cell, sArray);
    }
    /// ditto
    NeighList!(LowP) createOnConfig(Cell!(LowP) cell, SegmentedArray!(Vector!(LowP,3)) sArray){
        return new SimpleNeighList!(LowP)(cell, sArray);
    }
    void desc(CharSink s){
        s("<SimpleNeighListGen>");
    }
    mixin(serializeSome("dchem.SimpleNeighListGen","a neighboring list generator",""));
}
/// a neighboring list using the ellipsoid method (expensive unless cell<cutoff)
class SimpleNeighList(T):NeighList!(T){
    Cell!(T) cell;
    SegmentedArray!(Vector!(T,3)) sArray;

    KindRange kRange(){
        return sArray.kRange;
    }
    // default constructor just for serialization
    this(){}
    this(Cell!(T) cell, SegmentedArray!(Vector!(T,3))sArray){
        this.cell=cell;
        this.sArray=sArray;
    }
    /// loops on the neighbor of the point r with kind k
    int sloopOnNeighsRK(Vector!(T,3) r,KindIdx k,SafeT cutoff2,LoopPairs loopBody,
        size_t optSize=defaultOptimalBlockSize)
    {
        PIndex invalidP;
        auto arr2=sArray[KindRange(k,cast(KindIdx)(k+1))];
        foreach(ref i2,ref p2,ref l2,ref v;arr2.sLoop){
            int res=loopEllipsoidBounds!(T)(cell,r,v,cast(T)cutoff2,invalidP,0,p2,i2,loopBody);
            if (res!=0) return res;
        }
    }
    /// loops on the neighbor of p with kind k
    int sloopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2,LoopPairs loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        Vector!(T,3) v0=sArray[p1,i1];
        if (isNaN(cutoff2) || (!(k in sArray.kRange)) || sArray[k].length==0) return 0;
        SegmentedArray!(Vector!(T,3)) arr2;
        if (avoidDoubleCount){
            arr2=sArray[KindRange(k,cast(KindIdx)(k+1))];
            foreach(ref i2,ref p2,ref l2,ref v;arr2.sLoop){
                if (((p1.data+p2.data)&1UL)==cast(ulong)(p1.data>p2.data)) {
                    int res=loopEllipsoidBounds!(T)(cell,v0,v,cast(T)cutoff2,p1,i1,p2,i2,loopBody);
                    if (res!=0) return res;
                }
            }
        } else {
            foreach(ref i2,ref p2,ref l2,ref v;arr2.sLoop){
                int res=loopEllipsoidBounds!(T)(cell,v0,v,cast(T)cutoff2,p1,i1,p2,i2,loopBody);
                if (res!=0) return res;
            }
        }
        return 0;
    }
    /// loops on the neighbors between particles of kind k1 and k2
    int sloopOnNeighsKK(KindIdx k1,KindIdx k2,SafeT cutoff2,LoopPairs loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        if (isNaN(cutoff2) || sArray[k1].length==0 ||
            sArray[k2].length==0) return 0;
        if (avoidDoubleCount) {
            foreach(ref i1,ref p1,ref l1,ref v0;sArray[KindRange(k1,cast(KindIdx)(k1+1))].sLoop){
                foreach(ref i2,ref p2,ref l2,ref v;sArray[KindRange(k2,cast(KindIdx)(k2+1))].sLoop){
                    if (((p1.data+p2.data)&1UL)==cast(ulong)(p1.data>p2.data)){
                        int res=loopEllipsoidBounds!(T)(cell,v0,v,cast(T)cutoff2,p1,i1,p2,i2,loopBody);
                        if (res!=0) return res;
                    }
                }
            }
        } else {
            foreach(ref i1,ref p1,ref l1,ref v0;sArray[KindRange(k1,cast(KindIdx)(k1+1))].sLoop){
                foreach(ref i2,ref p2,ref l2,ref v;sArray[KindRange(k2,cast(KindIdx)(k2+1))].sLoop){
                    int res=loopEllipsoidBounds!(T)(cell,v0,v,cast(T)cutoff2,p1,i1,p2,i2,loopBody);
                    if (res!=0) return res;
                }
            }
        }
        return 0;
    }
    
    /// parallel loop on the neighbor of the point r with kind k
    int ploopOnNeighsRK(Vector!(T,3) r,KindIdx k,SafeT cutoff2,LoopPairs loopBody,
        size_t optSize=defaultOptimalBlockSize)
    {
        return sloopOnNeighsRK(r,k,cutoff2,loopBody,optSize);
    }
    /// parallel loop on the neighbor of p with kind k
    int ploopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2,LoopPairs loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return sloopOnNeighsPK(p1,i1,k,cutoff2,loopBody,avoidDoubleCount,optSize);
    }
    /// parallel loop on the neighbors between particles of kind k1 and k2
    int ploopOnNeighsKK(KindIdx k1,KindIdx k2,SafeT cutoff2, LoopPairs loopBody, bool avoidDoubleCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        return sloopOnNeighsKK(k1,k2,cutoff2,loopBody,avoidDoubleCount,optSize);
    }
    
    void release0(){
        // cleanup
    }
    
    mixin(serializeSome("dchem.SimpleNeighList!("~T.stringof~")","a neighboring list","cell|sArray"));
    mixin printOut!();
    mixin RefCountMixin!();
}
