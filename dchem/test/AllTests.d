module dchem.test.AllTests;
import blip.rtest.RTest;
import dchem.test.UtilTests: utilTests;
import dchem.test.SysTests: sysTests;
import dchem.test.OmgTests: omgTests;

/// collection of all the tests of dchem
TestCollection allTests()(TestCollection superColl){
    TestCollection coll=new TestCollection("dchem",
        __LINE__,__FILE__,superColl);
    utilTests(coll);
    omgTests(coll);
    sysTests!()(coll);
    return coll;
}
