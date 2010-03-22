/// the main dchem program
module Sampler;
//import dchem.input.RootInput;
import blip.io.Console;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import blip.t.core.stacktrace.TraceExceptions;
import blip.util.TraceAll;
import blip.text.TextParser;
import dchem.sys.SegmentedArray;
import dchem.input.RootInput;
import dchem.sampler.SinglePoint;
import dchem.calculator.AllCalculators;
import dchem.sampler.AllSamplers;

int main(char[][]args){
    if (args.length!=2){
        sout("Expected a single argument (inputfile)\n");
        return 1;
    }
    auto parser=new TextParser!(char)(toReaderT!(char)((new DataFileInput(args[1])).input));
    SegmentedArray!(real) sArr;
    auto rFile=new RootInput();
    rFile.readInput(parser,serr.call);
    auto mainInp="main" in rFile.knownNames;
    if (mainInp is null){
        sout("did not find 'main' in the input fields, nothing to do, stopping\n");
    } else {
        auto mainField=*mainInp;
        if (mainField.inputField !is null) mainField=mainField.inputField;
        auto samp=mainField.sampler;
        if (samp is null){
            sout("Error: 'main' should be a sampler, stopping\n");
        } else {
            samp.run();
        }
    }
    return 0;
}
