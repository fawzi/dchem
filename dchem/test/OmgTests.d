module dchem.test.OmgTests;
import xf.omg.core.LinearAlgebra;
import dchem.Common;
import blip.serialization.Serialization;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;

void m2NaTest(T)(SizedRandomNArray!(T,9)matBase){
    Matrix!(T,3,3) m1;
    
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            m1[i,j]=matBase.arr[3*i+j];
        }
    }
    auto m2=m2NA(&m1);
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            assert(m1[i,j]==m2[i,j]);
        }
    }
    auto m3=m2NAC(m1);
    for (int i=0;i<3;++i){
        for (int j=0;j<3;++j){
            assert(m1[i,j]==m3[i,j]);
        }
    }
    m2[1,1]=m2[1,1]+cscalar!(T,1);
    assert(m2[1,1]==m1[1,1]);
    assert(m3[1,1]!=m1[1,1]);
}

void v2NaTest(T)(SizedRandomNArray!(T,3)matBase){
    Vector!(T,3) v1;
    
    for (int i=0;i<3;++i){
        v1[i]=matBase.arr[i];
    }
    auto v2=v2NA(&v1);
    for (int i=0;i<3;++i){
        assert(v1[i]==v2[i]);
    }
    auto v3=v2NAC(v1);
    for (int i=0;i<3;++i){
        assert(v1[i]==v3[i]);
    }
    v2[1]=v2[1]+cscalar!(T,1);
    assert(v2[1]==v1[1]);
    assert(v3[1]!=v1[1]);
}

void addNaTests(T)(TestCollection coll){
    autoInitTst.testNoFailF("v2NaTestF!("~T.stringof~")",&v2NaTest!(T),
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("m2NaTest!("~T.stringof~")",&m2NaTest!(T),
        __LINE__,__FILE__,coll);
}

TestCollection omgTests(TestCollection superColl){
    TestCollection coll=new TestCollection("omg",__LINE__,__FILE__,superColl);
    addNaTests!(double)(coll);
    return coll;
}
    
