/// testing support to generate random cells
module dchem.test.sys.CellSupport;
import dchem.Common;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.narray.TestSupport;
import tango.math.Math;
import dchem.util.Rotate;
import dchem.sys.Cell;
import blip.io.Console;
import blip.io.BasicIO;
import blip.container.GrowableArray;

/// random orthorombic cell
class RandomOrthoCell{
    Cell cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        auto h=zerosR([3,3]);
        h[0,0]=r.uniformR2!(Real)(1.5,15.0);
        h[1,1]=r.uniformR2!(Real)(1.5,15.0);
        h[2,2]=r.uniformR2!(Real)(1.5,15.0);
        auto x0=emptyR(3);
        randomizeNArray(r.normalD(cast(Real)2.0),x0);
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)?1:0);
        }
        auto res=new typeof(this)();
        res.cell=new Cell(h,periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}
  
/// random cell with a along x axis, and b in the xy axis (h is upper triangular)
class RandomNromCell{
    Cell cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        auto param=emptyR(6);
        randomizeNArray(r.uniformD!(Real)(),param);
        param[Range(0,3)]*=13.5;
        param[Range(0,3)]+=1.5;
        param[Range(3,6)]*=150.0;
        param[Range(3,6)]+=15.0;
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        auto x0=emptyR(3);
        randomizeNArray(r.normalD(cast(Real)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell(cellParam2h(param),periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}
/// random cell
class RandomCell{
    Cell cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        scope param=emptyR([6]);
        randomizeNArray(r.uniformD!(Real)(),param);
        serr(collectAppender(delegate void(CharSink s){ s("param pre:");writeOut(s,param); s("\n");}));
        param[Range(0,3)]*=13.5;
        param[Range(0,3)]+=1.5;
        param[Range(3,6)]*=150.0;
        param[Range(3,6)]+=15.0;
        serr(collectAppender(delegate void(CharSink s){ s("param:");writeOut(s,param); s("\n");}));
        NArray!(Real,1) dir=emptyR([3]);
        randomizeNArray(r.normalD(cast(Real)1.0),dir);
        dir/=norm2(dir);
        scope h=cellParam2h(param);
        serr(collectAppender(delegate void(CharSink s){ s("h pre:");writeOut(s,h); s("\n");}));
        rotateVEi(dir,0,h);
        serr(collectAppender(delegate void(CharSink s){ s("h:");writeOut(s,h); s("\n");}));
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        scope x0=emptyR(3);
        randomizeNArray(r.normalD(cast(Real)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell(h,periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}
