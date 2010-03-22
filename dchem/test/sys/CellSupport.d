/// testing support to generate random cells
module dchem.test.sys.CellSupport;
import dchem.Common;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;
import tango.math.Math;
import dchem.util.Rotate;
import dchem.sys.Cell;
import blip.io.Console;
import blip.io.BasicIO;
import blip.container.GrowableArray;

/// random orthorombic cell
class RandomOrthoCell(T){
    Cell!(T) cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        auto h=zerosR([3,3]);
        h[0,0]=r.uniformR2!(T)(1.5,15.0);
        h[1,1]=r.uniformR2!(T)(1.5,15.0);
        h[2,2]=r.uniformR2!(T)(1.5,15.0);
        auto x0=emptyR(3);
        randomizeNArray(r.normalD(cast(T)2),x0);
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)?1:0);
        }
        auto res=new typeof(this)();
        res.cell=new Cell!(T)(h,periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}
  
/// random cell with a along x axis, and b in the xy axis (h is upper triangular)
class RandomNormCell(T){
    Cell!(T) cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        auto param=emptyR(6);
        randomizeNArray(r.uniformD!(T)(),param);
        param[Range(0,3)]*=13.5;
        param[Range(0,3)]+=1.5;
        param[Range(3,6)]*=150.0;
        param[Range(3,6)]+=15.0;
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        auto x0=emptyR(3);
        randomizeNArray(r.normalD(cast(T)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell(cellParamNArr2h(param),periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}
/// random cell
class RandomCell(T){
    Cell!(T) cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        scope param=emptyR([6]);
        randomizeNArray(r.uniformD!(T)(),param);
        serr(collectAppender(delegate void(CharSink s){ s("param pre:");writeOut(s,param); s("\n");}));
        param[Range(0,3)]*=13.5;
        param[Range(0,3)]+=1.5;
        param[Range(3,6)]*=150.0;
        param[Range(3,6)]+=15.0;
        serr(collectAppender(delegate void(CharSink s){ s("param:");writeOut(s,param); s("\n");}));
        NArray!(T,1) dir=emptyR([3]);
        randomizeNArray(r.normalD(cast(Real)1.0),dir);
        dir/=norm2(dir);
        scope h=cellParamNArr2h(param);
        scope hNA=m2NA(&h);
        serr(collectAppender(delegate void(CharSink s){ s("h pre:");writeOut(s,h); s("\n");}));
        rotateVEi(dir,0,hNA);
        serr(collectAppender(delegate void(CharSink s){ s("h:");writeOut(s,h); s("\n");}));
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        scope x0=emptyR(3);
        randomizeNArray(r.normalD(cast(T)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell!(T)(h,periodic,x0);
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}

/// test cell transfer functions
void testCell2Param(T)(T a,T b,T c,T alpha,T beta,T gamma,SizedRandomNArray!(T,3)dirA){
    auto h=cellParam2h(a,b,c,alpha,beta,gamma);
    auto param=h2CellParam(h);
    auto error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        throw new Exception(collectAppender(delegate void(CharSink s){
            s("Error too big:"); writeOut(s,error);
        }),__FILE__,__LINE__);
    }
    auto dir=dirA.arr;
    dir/=norm2(dir);
    scope hNA=m2NA(&h);
    auto h2=rotateVEi(dir,0,hNA);
    param=h2CellParam(h2);
    error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        throw new Exception(collectAppender(delegate void(CharSink s){
            s("Error2 too big:"); writeOut(s,error);
        }),__FILE__,__LINE__);
    }
}

void cell2ParamTests(T)(TestCollection coll){
  autoInitTst.testNoFail("testCell2Param!("~T.stringof~")",
      (T a,T b,T c,T alpha,T beta,T gamma,SizedRandomNArray!(T,3)dirA){ testCell2Param!(T)(a,b,c,alpha,beta,gamma,dirA);},
    __LINE__,__FILE__,coll);
}

/// collection of all the tests on the cell module
TestCollection cellTests(TestCollection superColl){
    TestCollection coll=new TestCollection("cell",
        __LINE__,__FILE__,superColl);
    cell2ParamTests!(double)(coll);
    return coll;
}

