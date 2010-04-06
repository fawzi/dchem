module dchem.test.sys.PIndexesTests;
import dchem.sys.PIndexes;
import blip.rtest.RTest;
import blip.io.BasicIO;
import blip.container.GrowableArray;

/// a random range (the range itself has no randomGenerate because normally one 
/// wants a range that is valid within a structure)
struct RandomKindRange{
    KindRange kRange;
    static RandomKindRange randomGenerate(Rand r){
        uint k1,k2;
        r(k1)(k2);
        k1=k1%30;
        k2=k2%30;
        if (k2<k1){
            if (k2==0){
                k2=k1;
            } else {
                auto t=k2;
                k2=k1;
                k1=t;
            }
        }
        RandomKindRange res;
        res.kRange=KindRange(cast(KindIdx)k1,cast(KindIdx)k2);
        return res;
    }
    void desc(CharSink s){
        writeOut(s,kRange);
    }
}

/// tests in and intersections of KindRanges
void testKindRangeIn(RandomKindRange rr1,RandomKindRange rr2,uint kk){
    KindRange r1=rr1.kRange, r2=rr2.kRange;
    auto r3=r1.intersect(r2);
    auto r4=r2.intersect(r1);
    assert(r3==r4,"intersect commutativity");
    KindIdx k=cast(KindIdx)(kk%30);
    if (r1 in r2){
        assert(r3==r1,"intersect equality");
        if (k in r1){
            assert(k in r2,"inclusion of sets => elements");
        }
    } else if (r2 in r1){
        assert(r3==r2,"inclusion equality");
        if (k in r2){
            assert(k in r3,"inclusion of sets => elements");
        }
    } else {
        if (k in r3){
            assert(k in r1 && k in r2,"inclusion of sets => elements");
        }
    }
    assert(!r1.isDummy,"dummy");
    assert(!r2.isDummy,"dummy");
    assert(!r3.isDummy,"dummy");
    assert(!r4.isDummy,"dummy");
    assert(r1.valid,"valid");
    assert(r2.valid,"valid");
    KindRange dummyR;
    assert(dummyR.isDummy,"dummyR");
}

/// a random PIndex (the PIndex itself has no randomGenerate because normally one 
/// wants a PIndex that is valid within a structure)
struct RandomPIndex(T){
    T pIndex;
    static RandomPIndex randomGenerate(Rand r){
        uint k,p;
        r(k)(p);
        k=k%30;
        p=p%100;
        RandomPIndex res;
        res.pIndex=T(cast(KindIdx)k,cast(ParticleIdx)p);
        return res;
    }
}

void testPIndexCreate(T)(uint kk,uint pp){
    auto k=kk%30;
    auto p=pp%100;
    auto v1=T(k,p);
    auto v2=T(cast(KindIdx)k,cast(ParticleIdx)p);
    auto v3=T(cast(short)k,cast(long)p);
    auto v4=T(cast(int)k,cast(int)p);
    auto v5=T(cast(ulong)k,cast(ulong)p);
    assert(cast(uint)v1.kind==k);
    assert(cast(uint)v1.particle==p);
    assert(cast(uint)v2.kind==k);
    assert(cast(uint)v2.particle==p);
    assert(v1==v2);
    assert(v1==v3);
    assert(v1==v4);
    assert(v1==v5);
    assert(v1.data==v2.data);
    assert(v1.data==v3.data);
    assert(v1.data==v4.data);
    assert(v1.data==v5.data);
    assert(v1.valid);
    T dummyP;
    assert(!dummyP.valid);
}

void testPIndexConvert(T,U)(RandomPIndex!(T) rI){
    auto v1=rI.pIndex;
    auto v1T=U(v1);
    assert(v1T.kind==v1.kind && v1T.particle==v1.particle);
    assert(v1T.data==v1.data);
    assert(T(v1T)==v1);
}

void testPIndexAdd(T)(RandomPIndex!(T) rI){
    auto pI=rI.pIndex;
    auto next1=pI;
    next1+=1;
    auto next2=pI;
    next2+=cast(ParticleIdx)1;
    auto next3=pI;
    next3+=cast(ulong)1;
    assert(next1.kind==pI.kind && next1.particle==cast(ParticleIdx)(1+pI.particle));
    assert(next1==next2);
    assert(next1==next3);
    assert(next1 != pI);
    next1+=-1;
    assert(next1==pI);
    next1+=cast(KindIdx)1;
    assert(next1.kind==cast(KindIdx)(1+pI.kind) && next1.particle==pI.particle);
}

void testPIndexSet(T)(RandomPIndex!(T) rI,uint k1,uint p1){
    k1=k1%30;
    p1=p1%100;
    auto i1=rI.pIndex;
    auto i2=i1;
    i2.kind=cast(KindIdx)k1;
    assert(i2.kind==cast(KindIdx)k1 && i2.particle==i1.particle);
    i2.kind=cast(int)k1+1;
    assert(i2.kind==cast(KindIdx)(k1+1) && i2.particle==i1.particle);
    i2.kind=i1.kind;
    assert(i1==i2);
    i2.particle=cast(ParticleIdx)p1;
    assert(i2.particle==cast(ParticleIdx)p1 && i1.kind==i2.kind);
    i2.particle=cast(int)p1+1;
    assert(i2.particle==cast(ParticleIdx)(p1+1) && i1.kind==i2.kind);
    i2.particle=i1.particle;
    assert(i1==i2);
}

TestCollection kindRangeTests(TestCollection superColl=null){
    TestCollection coll=new TestCollection("kindRange",
        __LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("in",
        &testKindRangeIn,__LINE__,__FILE__,coll);
    return coll;
}

TestCollection pIndexTests(T)(TestCollection superColl=null){
    TestCollection coll=new TestCollection(T.stringof,
        __LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("create",
        &testPIndexCreate!(T),__LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("convertPIndex",
        &testPIndexConvert!(T,PIndex),__LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("convertLocalPIndex",
        &testPIndexConvert!(T,LocalPIndex),__LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("add",
        &testPIndexAdd!(T),__LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("set",
        &testPIndexSet!(T),__LINE__,__FILE__,coll);
    return coll;
}

/// collection of all the tests on the cell module
TestCollection pIndexesTests()(TestCollection superColl=null){
    TestCollection coll=new TestCollection("PIndexes",
        __LINE__,__FILE__,superColl);
    kindRangeTests(coll);
    pIndexTests!(PIndex)(coll);
    pIndexTests!(LocalPIndex)(coll);
    return coll;
}

