module RandomTests;
import blip.narray.NArray;
import blip.rtest.RTest;
import dchem.test.AllTests: allTests;
import blip.math.random.Random;
import blip.io.Console;
version(NoTrace){} else { import blip.core.stacktrace.TraceExceptions; import blip.util.TraceAll; }

void main(char[][] args){
    serr(rand.toString());
    serr("\n");
    mainTestFun(args,allTests!()());
}
