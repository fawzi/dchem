/// context of the current process
module dchem.calculator.ProcContext;
import tango.math.random.Random;
import tango.io.vfs.FileFolder;
import tango.sys.Environment;
import blip.sync.UniqueNumber;

/// context of the current process, use the instance variable
class ProcContext{
    ulong id; /// unique id for this process
    VfsFolder baseDirectory; /// base directory where to write things
    UniqueNumber!(ulong) localId; /// a locally valid unique number
    
    this(){
        rand(id);
        localId=UniqueNumber!(ulong)(1);
        baseDirectory=new FileFolder(Environment.cwd());
    }
    /// global instance
    static ProcContext instance;
    static this(){
        instance=new ProcContext();
    }
}