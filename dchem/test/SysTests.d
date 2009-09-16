module dchem.test.SysTests;
import blip.rtest.RTest;
/+
/// test cell transfer functions
void testCell2Param(Real a,Real b,Real c,Real alpha,Real beta,Real gamma,SizedRandomNArray!(Real,3)dirA){
    auto h=cellParam2h(a,b,c,alpha,beta,gamma);
    auto param=h2CellParam(h);
    auto error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        Trace.formatln("error:{}",error);
        throw new Exception("Error too big",__FILE__,__LINE__);
    }
    auto dir=dirA.arr;
    dir/=norm2(dir);
    auto h2=rotateVEi(dir,0,h);
    param=h2CellParam(h2);
    error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        Trace.formatln("error2:{}",error);
        throw new Exception("Error2 too big",__FILE__,__LINE__);
    }
}
+/

/// collection of all the tests on the sys modules
TestCollection sysTests(TestCollection superColl){
    TestCollection coll=new TestCollection("sys",
        __LINE__,__FILE__,superColl);
    return coll;
}
