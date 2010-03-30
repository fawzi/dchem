module dchem.test.SysTests;
import blip.rtest.RTest;
import dchem.test.sys.CellSupport: cellTests;
import dchem.test.sys.SubMappingSupport: submapTests;
import dchem.test.sys.SegmentedArraySupport: segmentedArrayTests;


/// collection of all the tests on the sys modules
TestCollection sysTests()(TestCollection superColl){
    TestCollection coll=new TestCollection("sys",
        __LINE__,__FILE__,superColl);
    cellTests(coll);
    submapTests(coll);
    segmentedArrayTests!()(coll);
    return coll;
}
