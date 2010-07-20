module dchem.test.sampler.MinEExplorerTests;
import dchem.Common;
import blip.rtest.RTest;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import dchem.sampler.MinEExplorer;
import blip.math.Math: max;
import dchem.pnet.DirArray;

void testContains00(size_t v){
    size_t v2=v;
    for(int i=size_t.sizeof*4;i!=0;--i){
        if ((v2&0b11)==0){
            assert(FlagsArray.contains00(v)!=0,collectAppender(delegate void(CharSink s){
                s("error value ");s(printBits(v)); s(" contains 0, but contains00 did not catch it");
            }));
            return;
        }
        v2>>=2;
    }
    assert(FlagsArray.contains00(v)==0,collectAppender(delegate void(CharSink s){
        s("error value ");s(printBits(v)); s(" does not contain 0, but contains00 is not 0");
    }));
}

void testFlagsArrayFindFree(FlagsArray d,size_t s){
    size_t start=s;
    size_t refV=d.length;
    if (d.length==0){
        start=0;
        refV=0;
    } else {
        start=start%d.length;
        size_t endL=d.length-1;
        if (endL==start) endL=d.length;
        for (size_t i=0;i<d.length;++i){
            if (d[(i+start)%endL]==0){
                refV=(i+start)%endL;
                break;
            }
        }
    }
    if (refV==d.length && d.length>0 && d[d.length-1]==0) refV=d.length-1;
    auto findF=d.findFreeAndSet(start);
    assert(findF==refV,collectAppender(delegate void(CharSink s){
        dumper(s)("findFree gave a different result than expected:")(findF)(" vs ")(refV)
            (" for ")(d)(" starting at ")(start);
    }));
}

TestCollection flagsArrayTests(TestCollection superColl){
    TestCollection coll=new TestCollection("FlagsArray",__LINE__,__FILE__,superColl);
    autoInitTst.testNoFailF("contains00",&testContains00,
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("findFree",&testFlagsArrayFindFree,
        __LINE__,__FILE__,coll);
    return coll;
}

TestCollection minEETests()(TestCollection superColl){
    TestCollection coll=new TestCollection("minEE",__LINE__,__FILE__,superColl);
    flagsArrayTests!(double)(coll);
    return coll;
}
