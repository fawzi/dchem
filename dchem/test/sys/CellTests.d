/// testing support to generate random cells
module dchem.test.sys.CellTests;
import dchem.Common;
import blip.narray.NArray;
import blip.rtest.RTest;
import blip.test.narray.NArraySupport;
import tango.math.Math;
import dchem.util.Rotate;
import dchem.sys.Cell;
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

/// randome parameters for a cell (all angles in 15..165, sizes 1.5..13.5)
struct RandomCellParams(T){
    NArray!(T) param;
    static RandomCellParams randomGenerate(Rand r,ref bool acceptable){
        RandomCellParams res;
        auto param=emptyR([6]);
        randomizeNArray(r.uniformD!(T)(),param);
        param[Range(0,3)]*=13.5;
        param[Range(0,3)]+=1.5;
        param[Range(3,5)]*=150.0;
        param[Range(3,5)]+=15.0;
        auto minA=abs(param[3]-param[4])+14.99;
        auto maxA=min(165.0,min(360.0-15.0-param[3]+param[4],param[3]+param[4]-14.99));
        param[5]=cast(T)(minA+(maxA-minA)*param[5]);
        res.param=param;
        return res;
    }
    void desc(CharSink s){
        if (param is null){
            s("*null*");
        } else {
            param.printData(s);
        }
    }
}

/// random cell with a along x axis, and b in the xy axis (h is upper triangular)
class RandomNormCell(T){
    Cell!(T) cell;
    this(){}
    static typeof(this) randomGenerate(Rand r){
        RandomCellParams!(T) rParam;
        simpleRandom(rParam);
        param=rParam.param;
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
          periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        auto x0=emptyR(3);
        randomizeNArray(r.normalD(cast(T)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell(cellParamNArr2h(param),periodic,na2v3(x0));
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
    static typeof(this) randomGenerate(Rand r,ref bool acceptable){
        RandomCellParams!(T) rParam;
        acceptable = acceptable&&simpleRandom(r,rParam);
        scope param=rParam.param;
        NArray!(T,1) dir=emptyR([3]);
        randomizeNArray(r.normalD(cast(Real)1.0),dir);
        dir/=norm2(dir);
        auto h=cellParamNArr2h(param);
        scope hNA=m2NAC(h);
        auto h2=rotateVEi(dir,0,hNA);
        int[3] periodic;
        for (int idim=0;idim<3;++idim){
            periodic[idim]=(r.uniform!(bool)()?1:0);
        }
        scope x0=emptyR(3);
        randomizeNArray(r.normalD(cast(T)2.0),x0);
        auto res=new typeof(this)();
        res.cell=new Cell!(T)(h,periodic,na2v3(x0));
        return res;
    }
    void desc(CharSink s){
        cell.desc(s);
    }
}

/// test cell transfer functions
void testCell2Param(T)(RandomCellParams!(T)rParams,SizedRandomNArray!(T,3)dirA){
    T angErr(real a1){
        auto redE=a1/360;
        return abs(redE-floor(redE+0.5)); // does not multiply by 360!
    }
    auto p=rParams.param;
    auto h=cellParamNArr2h(rParams.param);
    auto param=h2CellParam(h);
    auto error=abs(param[0]-p[0])+abs(param[1]-p[1])+abs(param[2]-p[2])+angErr(param[3]-p[3])
        +angErr(param[4]-p[4])+angErr(param[5]-p[5]);
    if (error>1.e-6) {
        throw new Exception(collectAppender(delegate void(CharSink s){
            s("Error too big:"); writeOut(s,error);
            s("calcParam:"); writeOut(s,param);
            s("h:"); writeOut(s,h);
        }),__FILE__,__LINE__);
    }
    auto dir=dirA.arr;
    dir/=norm2(dir);
    scope hNA=m2NAC(h);
    auto h2=rotateVEi(dir,0,hNA);
    param=h2CellParam(h2);
    error=abs(param[0]-p[0])+abs(param[1]-p[1])+abs(param[2]-p[2])+angErr(param[3]-p[3])
        +angErr(param[4]-p[4])+angErr(param[5]-p[5]);
    if (error>1.e-6) {
        throw new Exception(collectAppender(delegate void(CharSink s){
            s("Error2 too big:"); writeOut(s,error);
            s("calcParam:"); writeOut(s,param);
        }),__FILE__,__LINE__);
    }
}

void cell2ParamTests(T)(TestCollection coll){
  autoInitTst.testNoFailF("testCell2Param!("~T.stringof~")",
      &testCell2Param!(T),__LINE__,__FILE__,coll);
}

/// collection of all the tests on the cell module
TestCollection cellTests(TestCollection superColl){
    TestCollection coll=new TestCollection("cell",
        __LINE__,__FILE__,superColl);
    cell2ParamTests!(double)(coll);
    return coll;
}

