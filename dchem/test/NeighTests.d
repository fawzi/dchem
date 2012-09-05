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

class NeighListTester(T) {
    int[NeighId] neighList;
    Cell!(T) cell;
    int expectedValue;
    SegmentedArray!(Vector!(T,3)) pos;
    SimpleNeighList!(T) sNeighs;
    HierarchicalSort!(T) hNeighs;
    PIndex iPart1;
    uint iPos;
    PIndex p1;
    KindIdx kind1, kind2;
    RandomSegmentedArray!(T) rArray;
    SafeT cutoff2;
    bool useLock=false;
    bool weakVerify=false;
    char[] context;
    
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
            if (weakVerify && ret==0) {
                if (useLock){
                    synchronized(this){
                        neighList[nId] += expectedValue-1;
                    }
                } else {
                    neighList[nId] += expectedValue-1;
                }
            } else {
                throw new Exception(collectAppender(delegate void(CharSink s){
                    dumper(s)("doLoop: Error for neighs ")(neighs)(" expected ")(expectedValue)
                        (", had ")(ret)(" in ")(context);
                }),__FILE__,__LINE__);
            }
        }
        return 0;
    }
    
    void verifyValues(){
        foreach (k,v; neighs){
            if (v!=expectedValue) {
                if (weakVerify && v==expectedValue-1) {
                    neighList.remove[k];
                } else {
                    throw new Exception(collectAppender(delegate void(CharSink s){
                        dumper(s)("verifyValues: Error for neighs ")(k)(" expected ")(expectedValue)
                            (", had ")(v)(" in ")(context);
                    }),__FILE__,__LINE__);
                }
            }
        }
    }
    
    this(RandomCell!(T) rcell, RandomSegmentedArray!(T) rArray) {
        this.cell=rcell.cell;
        this.pos=rArray.randomArray();
        auto sGenTest=new SimpleNeighGen();
        sNeighs=sGenTest.createOnConfig(cell,pos);
        auto hGenTest=new HierarchicalSortGen();
        hNeighs=hGenTest.createOnConfig(cell,pos);
    }
    
    void tests(uint iPart1, uint iPos1, uint ikind2, Vector!(T,3) randomR, SafeT cutoff2,
        SizeLikeNumber!(4,1,-3) bs, bool avoidDoubleCount){
        this.iPart1=iPart1;
        this.blockSize=bs.value;
        auto particles=this.pos.arrayStruct.submapping().sortedPIndex();
        p1=particles[iPos1%particles.length];
        this.iPos1=(iPos1%particles.length)%this.pos.arrayStruct.nParticles(p1.kind);
        this.avoidDoubleCount = avoidDoubleCount;
        this.kind1=p1.kind;
        this.kind2=cast(KindIdx)(particles.kRange.kStart+(ikind2%(particles.kRange.kEnd-particles.kRange.kStart)));
        this.cutoff2=cutoff2;
        
        useLock=false;
        context="sloopOnNeighsPK_1";
        sNeighs.sloopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        context="sloopOnNeighsPK_2";
        sNeighs.sloopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        useLock=true;
        context="ploopOnNeighsPK_1";
        hNeighs.ploopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,avoidDoubleCount, blockSize);
        ++expectedValue;
        verifyValues();
        context="ploopOnNeighsPK_2";
        hNeighs.ploopOnNeighsPK(p1, iPos1,kind2,cutoff2, &doLoop,avoidDoubleCount, blockSize);
        ++expectedValue;
        verifyValues();

        int[NeighId] neighListConf;
        foreach(n,i;neighList){
            NeighId nConv=n;
            PIndex dummy;
            nConv.p1=dummy;
            nConv.i1=0;
            neighListConf[nConv]=i;
        }
        neighList.clear();
        neighList=neighListConf;

        useLock=false;
        weakVerify=true;
        context="sloopOnNeighsRK_1";
        sNeighs.sloopOnNeighsRK(r1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        weakVerify=false;
        context="sloopOnNeighsRK_2";
        sNeighs.sloopOnNeighsRK(r1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        useLock=true;
        context="ploopOnNeighsRK_1";
        hNeighs.ploopOnNeighsRK(r1,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();
        context="ploopOnNeighsRK_2";
        hNeighs.ploopOnNeighsRK(r1,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();

        neighList.clear();
        useLock=false;
        context="sloopOnNeighsRK_1b";
        sNeighs.sloopOnNeighsRK(randomR,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        context="sloopOnNeighsRK_2b";
        sNeighs.sloopOnNeighsRK(randomR,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        useLock=true;
        context="ploopOnNeighsRK_1b";
        hNeighs.ploopOnNeighsRK(randomR,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();
        context="ploopOnNeighsRK_2b";
        hNeighs.ploopOnNeighsRK(randomR,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();

        neighList.clear();
        useLock=false;
        context="sloopOnNeighsKK_1";
        sNeighs.sloopOnNeighsKK(kind1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        context="sloopOnNeighsKK_2";
        sNeighs.sloopOnNeighsKK(kind1,kind2,cutoff2, &doLoop,avoidDoubleCount);
        ++expectedValue;
        verifyValues();
        useLock=true;
        context="ploopOnNeighsKK_1";
        hNeighs.ploopOnNeighsKK(kind1,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();
        context="ploopOnNeighsKK_2";
        hNeighs.ploopOnNeighsKK(kind1,kind2,cutoff2, &doLoop,avoidDoubleCount,blockSize);
        ++expectedValue;
        verifyValues();
        
        neighList.clear();
    }
}

void doNeigListTests(T)(RandomCell!(T) rcell, RandomSegmentedArray!(T) rArray,
    uint iPart1, uint iPos1, uint ikind2, Vector!(T,3) randomR, SafeT cutoff2,
    SizeLikeNumber!(4,1,-3) bs, bool avoidDoubleCount)
{
    auto tester=new NeighListTester!(T)();
    tester.tests()
}
