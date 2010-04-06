module dchem.test.sys.SysTests;
import blip.rtest.RTest;
import dchem.test.sys.CellTests: cellTests;
import dchem.test.sys.SubMappingTests: submapTests;
import dchem.test.sys.SegmentedArrayTests: segmentedArrayTests;
import dchem.test.sys.PIndexesTests: pIndexesTests;


/// collection of all the tests on the sys modules
TestCollection sysTests()(TestCollection superColl){
    TestCollection coll=new TestCollection("sys",
        __LINE__,__FILE__,superColl);
    cellTests(coll);
    pIndexesTests!()(coll);
    submapTests(coll);
    segmentedArrayTests!()(coll);
    return coll;
}
