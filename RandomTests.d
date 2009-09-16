module RandomTests;
import blip.narray.NArray;
import blip.rtest.RTest;
import dchem.test.AllTests: allTests;
import tango.math.random.Random;
import tango.util.log.Trace;

void main(char[][] args){
    Trace.formatln("{}",rand.toString());
    mainTestFun(args,allTests(null));
}
