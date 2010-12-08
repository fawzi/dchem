/// the main dchem program
module Sampler;
//import dchem.input.RootInput;
import blip.io.Console;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import blip.core.stacktrace.TraceExceptions;
import blip.util.TraceAll;
import blip.text.TextParser;
import dchem.sys.SegmentedArray;
import dchem.input.RootInput;
import dchem.sampler.SinglePoint;
import dchem.calculator.AllCalculators;
import dchem.sampler.AllSamplers;
import blip.io.EventWatcher;
import blip.parallel.rpc.Rpc;
import blip.parallel.mpi.Mpi;
import dchem.pnet.WorkAsker;
import blip.stdc.stdlib:exit;

int main(char[][]args){
    ProtocolHandler.defaultProtocol.startServer(false); // starting the rpc server...
    if (args.length!=2){
        sout("Expected a single argument (inputfile)\n");
        return 1;
    }
    auto parser=new TextParser!(char)(toReaderT!(char)((new DataFileInput(args[1])).file.input));
    SegmentedArray!(real) sArr;
    auto rFile=new RootInput();
    auto inputOk=rFile.readInput(parser,serr.call);
    if (!inputOk){
        sout("error while reading input!\n");
        exit(1);
    }
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
            samp.run(mpiWorld,sout.call);
        }
    }
    noToutWatcher.stopLoop();
    exit(0);
    return 0;
}
