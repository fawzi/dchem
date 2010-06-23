/// various utilities to execute external commands
module dchem.util.ExecCmd;
import tango.sys.Process;
import blip.io.IOArray;
import blip.io.NullStream;

version(Posix){
    import tango.stdc.posix.signal;

    static this(){
        sigaction_t sigpipe = void;
        (cast(byte*) &sigpipe)[0 .. sigaction_t.sizeof] = 0;
        sigpipe.sa_handler=SIG_IGN;
        auto status = sigaction( SIGPIPE, &sigpipe, null );
        if (status !=0){
            throw new Exception("could not ignore SIGPIPE",__FILE__,__LINE__);
        }
    }
}

/// creates a process that will be able to execute the given command
Process cmd(char[] cmd,char[] baseDir=""){
    auto p=new Process(cmd);
    if (baseDir.length){
        p.workDir=baseDir;
    }
    return p;
}
/// gets output (and errors) of the execution of the given process
char[] getOutput(Process p,out int status){
    auto buf=new IOArray(512,512);
    p.redirect=Redirect.Output|Redirect.ErrorToOutput;
    p.execute();
    buf.copy(p.stdout);
    auto res=p.wait();
    status=res.status;
    switch (res.reason)
    {
        case res.Exit:
            break;
        case res.Signal,res.Stop,res.Continue,res.Error:
        default:
            if (status==0) status=-1;
            break;
    }
    p.stdout.close();
    return cast(char[])buf.slice();
}
/// executes the given process, discards any output and return the status
int execProcess(Process p){
    auto buf=nullStream();
    p.redirect=Redirect.Output|Redirect.ErrorToOutput;
    p.execute();
    buf.copy(p.stdout);
    auto res=p.wait();
    auto status=res.status;
    switch (res.reason)
    {
        case res.Exit:
            break;
        case res.Signal,res.Stop,res.Continue,res.Error:
        default:
            if (status==0) status=-1;
            break;
    }
    p.stdout.close();
    return status;
}

