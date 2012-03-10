module dchem.test.NeighTests;
import dchem.Common;
import dchem.sys.Cell;
import dchem.neigh.SimpleNeigh;
import blip.rtest.RTest;
import dchem.test.sys.CellTests;
import dchem.sys.PIndexes;
import blip.util.Convert;
import dchem.neigh.NeighModels;
import blip.math.Math;
import dchem.neigh.HierarchicalSort;

void testEllipsoidBoundsT(T)(RandomCell!(T) cell, Vector!(T,3) r1, Vector!(T,3) r2) {
    loopEllipsoidBounds!(T,true)(cell.cell,r1,r2,convertTo!(T)(2),
    PIndex(0UL),0,PIndex(1UL),1,delegate int(ref NeighPair!(T)){ return 0; });
}

TestCollection neighTests(TestCollection superColl){
    TestCollection coll=new TestCollection("neighs",
        __LINE__,__FILE__,superColl);
    alias float T;
    autoInitTst.testNoFailF("testEllipsoidBoundsT!("~T.stringof~")",
        &testEllipsoidBoundsT!(T),__LINE__,__FILE__,coll);
    return coll;
}

struct NeighId{
    PIndex p1;
    index_type i1;
    PIndex p2;
    index_type i2;
    index_type i,j,k;
    
    equals_t opEquals(ref NeighId o){
        return p1==o.p1 && i1==o.i1 && p2==o.p2 && i2==o.i2
            && i==o.i && j==o.j && k==o.k;
    }
    
    int opCmp(ref NeighId o){
        if (p1.data < o.p1.data) return -1;
        if (p1.data > o.p1.data) return 1;
        if (i1 < o.i1) return -1;
        if (i1 > o.i1) return 1;
        if (p2.data < o.p2.data) return -1;
        if (p2.data > o.p2.data) return 1;
        if (i2 < o.i2) return -1;
        if (i2 > o.i2) return 1;
        if (i < o.i) return -1;
        if (i > o.i) return 1;
        if (j < o.j) return -1;
        if (j > o.j) return 1;
        if (k < o.k) return -1;
        if (k > o.k) return 1;
        return 0;
    }
    
    static NeighId opCall(T,U)(ref NeighPair!(T) neighs, Cell!(U) cell,real maxerr=0.01){
        NeighId res;
        res.p1=neighs.p1;
        res.i1=neighs.i1;
        res.p2=neighs.p2;
        res.i2=neighs.i2;
        auto iVect=cell.hInv*neighs.p12;
        res.i=floor(iVect.x+0.5);
        res.j=floor(iVect.y+0.5);
        res.k=floor(iVect.z+0.5);
        if (pow2(iVect.x-i)+pow2(iVect.y-j)+pow2(iVect.z-k)>maxerr*maxerr)
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("Lattice error too big for ")(neighs)(": ")(iVect)(" vs [")(i)(", ")(j)(", ")(k)("] error=")
                    (sqrt(pow2(iVect.x-i)+pow2(iVect.y-j)+pow2(iVect.z-k)))(">")(maxerr);
            }),__FILE__,__LINE__);
        return res;
    }
}

class NeighList(T) {
    int[NeighId] neighList;
    Cell!(T) cell;
    int expectedValue;
    SegmentedArray!(Vector!(T,3)) pos;
    SimpleNeighList!(T) sNeighs;
    HierarchicalSort!(T) hNeighs;
    PIndex iPart1;
    uint iPos;
    RandomSegmentedArray!(T) rArray;
    SafeT cutoff2;
    bool useLock=false;
    
    int addNeigh(ref NeighId nId){
        auto v = nId in neighList;
        int oldVal;
        if (v) {
            oldVal= *v;
        } else {
            oldVal=0;
        }
        neighList[nId]=oldVal+1;
        return oldVal;
    }
    
    int doLoop(U)(ref NeighPair!(U) neighs){
        auto r1=pos[neighs.p1,neighs.i1];
        auto r2=pos[neighs.p2,neighs.i2];
        NeighId nId=NeighId(neighs,cell,r1,r2);
        int ret;
        if (useLock){
            synchronized(this){
                ret=addNeigh(nId);
            }
        } else {
            ret=addNeigh(nId);
        }
        if (ret != expectedValue) {
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("doLoop: Error for neighs ")(neighs)(" expected ")(expectedValue)
                    (", had ")(ret);
            }),__FILE__,__LINE__);
        }
        return 0;
    }
    
    void verifyValues(){
        foreach (k,v; neighs){
            if (v!=expectedValue) {
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("verifyValues: Error for neighs ")(k)(" expected ")(expectedValue)
                        (", had ")(v);
                }),__FILE__,__LINE__);
            }
        }
    }
    
    void setup() {
        auto sGenTest=new SimpleNeighGen();
        sNeighs=sGenTest.createOnConfig(cell,pos);
        auto hGenTest=new HierarchicalSortGen();
        hNeighs=hGenTest.createOnConfig(cell,pos);
    }
    
    void tests(uint iPart1, uint iPos1, RandomSegmentedArray!(T) rArray, SafeT cutoff2){
        this.rArray=rArray;
        this.iPart1=iPart1;
        this.pos=rArray.randomArray();
        auto particless=this.pos.arrayStruct.submapping().sortedPIndex();
        p1=particles[iPos1%particles.length];
        this.iPos1=iPos1%this.pos.arrayStruct.nParticles(p1.kind);
        this.cutoff2=cutoff2;
        useLock=false;
        sNeighs.sloopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        sNeighs.sloopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        useLock=true;
        hNeighs.ploopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        hNeighs.ploopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        neighList.clear();
        sNeighs.sloopOnNeighsKK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        sNeighs.sloopOnNeighsKK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        hNeighs.sloopOnNeighsKK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        hNeighs.sloopOnNeighsKK(p1, iPos1,kind2,cutoff2, &doLoop,false);
        ++expectedValue;
        verifyValues();
        
        /// loops on the neighbor of the point r with kind k (p1 will be invalid, i1=0)
        int sloopOnNeighsRK(Vector!(T,3) r,KindIdx k,SafeT cutoff2,LoopPairs loopBody,
            size_t optSize=defaultOptimalBlockSize);
        /// loops on the neighbor of p with kind k
        int sloopOnNeighsPK(PIndex p1,index_type i1,KindIdx k,SafeT cutoff2,LoopPairs loopBody, bool avoidDoubleCount=false,
            size_t optSize=defaultOptimalBlockSize);
        /// loops on the neighbors between particles of kind k1 and k2
        int sloopOnNeighsKK(KindIdx k1,KindIdx k2,SafeT cutoff2,LoopPairs loopBody, bool avoidDoubleCount=false,
            size_t optSize=defaultOptimalBlockSize);
        
    }
}