module dchem.input.ReadCoordFile;
import dchem.input.RootInput;
import dchem.sys.ParticleSys;
import dchem.sys.SegmentedArray;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import dchem.input.ReadIn;
import dchem.Common;
import blip.util.TangoLog;
import blip.core.Variant;
import blip.io.StreamConverters;
import tango.io.stream.DataFile;
import tango.io.FilePath;
import dchem.input.ReadIn2PSys;
import dchem.sys.PIndexes;
import blip.core.Array;

/// configuration
interface Config:InputElement{
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem();
    /// the configuration as ParticleSys (might be cached)
    ParticleSys!(Real) particleSysReal();
    /// the configuration as low precision ParticleSys (might be cached)
    ParticleSys!(LowP) particleSysLowP();
    /// drops the cached readSystem/particleSys
    void clear();
}

/// reads an xyz file
class FileConfig:Config{
    InputField _myField;
    char[] fileName;
    long frame;
    char[] format;
    ReadSystem readSys;
    ParticleSys!(Real) pSysReal;
    ParticleSys!(LowP) pSysLowP;
    Real[][] cell;
    
    /// the read configuration as ReadSystem, this might trigger the real read, and might cache the result
    ReadSystem readSystem(){
        if (readSys is null){
            scope dfile=new DataFileInput(fileName);
            scope file=new MultiInput(dfile);
            readSys=readFrame(file,format,frame);
            file.shutdownInput();
            if (cell.length!=0){
                if (cell.length!=3) throw new Exception("invalid cell size (should be 3x3)",__FILE__,__LINE__);
                for (int i=0;i<3;++i) {
                    if (cell[i].length!=3) throw new Exception("invalid cell size (should be 3x3)",__FILE__,__LINE__);
                    for (int j=0;j<3;++j){
                        readSys.cell[i][j]=cell[i][j];
                    }
                }
            }
        }
        return readSys;
    }
    /// the configuration as ParticleSys (might be cached)
    ParticleSys!(Real) particleSysReal(){
        if (pSysReal is null){
            pSysReal=readIn2PSys!(Real)(readSystem());
        }
        return pSysReal.dup(PSDupLevel.All);
    }
    ParticleSys!(LowP) particleSysLowP(){
        if (pSysLowP is null){
            pSysLowP=readIn2PSys!(LowP)(readSystem());
        }
        return pSysLowP.dup(PSDupLevel.All);
    }
    /// drops the cached readSystem/particleSys
    void clear(){
        readSys=null;
        pSysReal=null;
    }
    mixin myFieldMixin!();
    bool verify(CharSink logger){
        bool res=true;
        auto w=dumper(logger);
        if (fileName.length==0){
            w("Error: fileName is empty in field ")(myFieldName)("\n");
            res=false;
        } else if (!FilePath(fileName).exists()){
            w("Warning: fileName points to a non existing file '")(fileName)("' in field ")(myFieldName)(", continuing expecting it to be created in future\n");
        };
        if (format.length==0){
            if (fileName.length>4){
                switch(fileName[$-4..$]){
                case ".xyz":
                    format="xyz";
                    break;
                case ".car":
                    format="car";
                    break;
                case ".pdb":
                    format="pdb";
                    break;
                default:
                    break;
                }
            }
            if (format.length==0) {
                w("Error: format not given, an could not derive it from fileName in field ")(myFieldName)("\n");
                res=false;
            }
        } else {
            if (find(readFrameFormats,format)==readFrameFormats.length){
                w("Error: unknown format '")(format)("' given in field ")(myFieldName)(", known formats are ");
                foreach(i,f;readFrameFormats){
                    if (i!=0) w(", ");
                    w(f);
                }
                w("\n");
                res=false;
            }
        }
        return res;
    }
    mixin(serializeSome("dchem.FileConfig",`a configuration that is read from a file`,
    `fileName: the filename where to reat the configuration
    format: the format of the file(xyz,car,pdb)
    frame: if a frame different from the first one should be read (-1 means the last one)
    cell: the cell matrix h (if present overrides the one that might be in the file)`));
    mixin printOut!();
}
