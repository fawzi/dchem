/// the main dchem program
module Sampler;
import blip.io.Console;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import blip.core.stacktrace.TraceExceptions;
import blip.util.TraceAll;
import blip.text.TextParser;
import dchem.sys.SegmentedArray;
import dchem.sys.AllConstraints;
import dchem.input.RootInput;
import dchem.sampler.SinglePoint;
import dchem.calculator.AllCalculators;
import dchem.sampler.AllSamplers;
import blip.io.EventWatcher;
import blip.parallel.rpc.Rpc;
import blip.parallel.mpi.Mpi;
import dchem.pnet.WorkAsker;
import blip.core.Thread;
import dchem.InputDesc;
import blip.stdc.stdlib:exit;
import blip.io.FileStream;

struct SamplerRunner{
    Sampler sampler;
    LinearComm paraEnv;
    CharSink log;
    
    void run(){
        sampler.run(paraEnv,log);
    }
}

int main(char[][]args){
    try{
        bool strict=false;
        foreach_reverse (i,a;args){
            bool drop=true;
            if (a.length>7 && a[0..7]=="--port="){
                auto stcp=cast(StcpProtocolHandler)ProtocolHandler.defaultProtocol;
                if (stcp!is null) stcp.port=a[7..$].dup;// non aligned... thus dup
            } else if (a=="--strict") {
                strict=true;
            } else if (a=="--help"){
                sout(args[0])(` [--help] [--port=portNr] [--strict] [--html-doc] inputFile

    --help: prints this help
    --port=portNr: starts the rpc server on portNr instead of 50000
    --strict: stops if the server could not be started
    --html-doc: generates a reference of the input in html format to the doc.html file`);
                sout("\n");
                exit(0);
            } else if (a=="--html-doc"){
                auto doc=new InputDesc();
                doc.loadTypes();
                auto outF=outfileStr("doc.html",WriteMode.WriteClear);
                try{
                    doc.htmlDesc(outF.charSink);
                } finally {
                    outF.close();
                }
                exit(0);
            } else {
                drop=false;
            }
            if (drop){
                for(size_t j=i+1;j<args.length;++j){ 
                    args[j-1]=args[j]; 
                }
                args=args[0..$-1];
            }
        }
        ProtocolHandler.defaultProtocol.startServer(strict); // starting the rpc server...
        sout("local rpc server started with url ")(ProtocolHandler.defaultProtocol.handlerUrl)("\n");
        if (args.length!=2){
            sout("Expected a single argument (inputfile)\n");
            exit(1);
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
                auto closure=new SamplerRunner;
                closure.sampler=samp;
                closure.paraEnv=mpiWorld;
                closure.log=sout.call;
                Task("mainSampler",&closure.run).autorelease.executeNow();
                sout("mainSampler finished!\n");
            }
            noToutWatcher.stopLoop();
        }
    } catch (Exception e){
        serr("Failure in main:")(e)("\n");
        Thread.sleep(0.3);
        exit(1);
    }
    exit(0);
    return 0;
}
