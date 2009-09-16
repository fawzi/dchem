module testUtilRotate;
import blip.narray.NArray;
import blip.rtest.RTest;
import dchem.util.Rotate;
import tango.math.random.Random;
import tango.util.log.Trace;
import blip.narray.TestSupport;

// tests rotation routine v1 v2 in 3D
void doTestRotV1V2(T,int ndim)(SizedRandomNArray!(T,ndim) v1,SizedRandomNArray!(T,ndim) v2){
  auto n1=v1.arr/norm2(v1.arr);
  auto n2=v2.arr/norm2(v2.arr);

  auto rotM=rotateVV!(NArray!(T,1),NArray!(T,2))(n1,n2,eye!(T)(ndim));
  auto err=minFeqrel2(dot(rotM,n1),n2);
  auto err2=minFeqrel2(dot(rotM,rotM.T),eye!(T)(ndim));
  if (err<T.mant_dig/4*3-10){
    Trace.formatln("err: {}/{}",err,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
  if (err2<T.mant_dig/4*3-10){
    Trace.formatln("err2: {}/{}",err2,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
}
/// tests rotation routines v1 j, i v2
void doTestRotV1J(T,int ndim)(SizedRandomNArray!(T,ndim) v1,size_t i){
  size_t idx=i%3;
  auto n1=v1.arr/norm2(v1.arr);
  auto n2=zeros!(T)(ndim);
  n2[idx]=cast(T)1;

  auto rotM_2=rotateEiV!(NArray!(T,1),NArray!(T,2))(idx,n1,eye!(T)(ndim));
  auto err_2=minFeqrel2(dot(rotM_2,n2),n1);
  auto err2_2=minFeqrel2(dot(rotM_2,rotM_2.T),eye!(T)(ndim));
  if (err_2<T.mant_dig/4*3-10){
    Trace.formatln("err: {}/{}",err_2,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
  if (err2_2<T.mant_dig/4*3-10){
    Trace.formatln("err2: {}/{}",err2_2,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }

  auto rotM=rotateVEi!(NArray!(T,1),NArray!(T,2))(n1,idx,eye!(T)(ndim));
  auto err=minFeqrel2(dot(rotM,n1),n2);
  auto err2=minFeqrel2(dot(rotM,rotM.T),eye!(T)(ndim));
  if (err<T.mant_dig/4*3-10){
    Trace.formatln("err: {}/{}",err,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
  if (err2<T.mant_dig/4*3-10){
    Trace.formatln("err2: {}/{}",err2,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
  
  auto err_3=minFeqrel2(rotM.T,rotM_2);
  if (err_3<T.mant_dig/4*3-10){
    Trace.formatln("err: {}/{}",err_3,T.mant_dig);
    throw new Exception("inconsistent V1Ei EiV1",__FILE__,__LINE__);
  }
}

void addRotTstToCollection(T,int ndim)(TestCollection coll){
  autoInitTst.testNoFail("rotV1V2!("~T.stringof~")",(SizedRandomNArray!(T,ndim)v1,SizedRandomNArray!(T,ndim)v2){doTestRotV1V2!(T)(v1,v2);},
    __LINE__,__FILE__,coll);
   autoInitTst.testNoFail("rotV1J!("~T.stringof~")",(SizedRandomNArray!(T,ndim)v1,size_t i){doTestRotV1J!(T)(v1,i);},
    __LINE__,__FILE__,coll);
}

TestCollection rotateTests(TestCollection superColl){
    TestCollection coll=new TestCollection("rotV1V2",
        __LINE__,__FILE__,superColl);
    addRotTstToCollection!(float,3)(coll);
    addRotTstToCollection!(double,3)(coll);
    addRotTstToCollection!(real,3)(coll);
    return coll;
}

void main(char[][] args){
    Stdout(rand.toString()).newline;
    mainTestFun(args,rotateTests(null));
}
