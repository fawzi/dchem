module dchem.test.UtilTests;
import blip.narray.NArray;
import blip.rtest.RTest;
import dchem.util.Rotate;
import tango.math.random.Random;
import tango.util.log.Trace;
import blip.test.narray.NArraySupport;

// tests rotation routine v1 v2 in 3D
void doTestRotV1V2(T,int ndim)(SizedRandomNArray!(T,ndim) v1,SizedRandomNArray!(T,ndim) v2){
  auto n1=v1.arr/norm2(v1.arr);
  auto n2=v2.arr/norm2(v2.arr);

  auto rotM=rotateVV!(NArray!(T,1),NArray!(T,2))(n1,n2,eye!(T)(ndim));
  auto err=minFeqrel2(dot(rotM,n1),n2);
  auto err2=minFeqrel2(dot(rotM,rotM.T),eye!(T)(ndim));
  if (err<T.mant_dig/4*3-11){
    Trace.formatln("err: {}/{}",err,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
  if (err2<T.mant_dig/4*3-11){
    Trace.formatln("err2: {}/{}",err2,T.mant_dig);
    throw new Exception("incorrect rotation",__FILE__,__LINE__);
  }
}
/// invariant v1v2
void invariantRotV1V2(T,int ndim,U=SizedRandomNArray!(T,ndim,ndim))(SizedRandomNArray!(T,ndim) v1,SizedRandomNArray!(T,ndim) v2,U m){
  auto n1=v1.arr/norm2(v1.arr);
  auto n2=v2.arr/norm2(v2.arr);

  auto rotM=rotateVV!(NArray!(T,1),NArray!(T,2))(n1,n2,m.arr.dup);
  auto rotM2=rotateVV!(NArray!(T,1),NArray!(T,2))(n2,n1,rotM);
  auto err=minFeqrel2(rotM2,m.arr);
  if (err<T.mant_dig/4*3-11){
    Trace.formatln("err: {}/{}",err,T.mant_dig);
    throw new Exception("non invariant",__FILE__,__LINE__);
  }
}

/// tests rotation routines v1 j, i v2
void doTestRotV1J(T,int ndim)(SizedRandomNArray!(T,ndim) v1,size_t i,T scale){
  size_t idx=i%ndim;
  auto n1=v1.arr/norm2(v1.arr);
  auto n2=zeros!(T)(ndim);
  n2[idx]=cast(T)1;

  {
      auto rotM_2=rotateEiV!(NArray!(T,1),NArray!(T,2))(idx,n1,eye!(T)(ndim));
      auto err_2=minFeqrel2(dot(rotM_2,n2),n1);
      auto err2_2=minFeqrel2(dot(rotM_2,rotM_2.T),eye!(T)(ndim));
      if (err_2<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err_2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
      if (err2_2<T.mant_dig/4*3-11){
        Trace.formatln("err2: {}/{}",err2_2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }

      auto rotM=rotateVEi!(NArray!(T,1),NArray!(T,2))(n1,idx,eye!(T)(ndim));
      auto err=minFeqrel2(dot(rotM,n1),n2);
      auto err2=minFeqrel2(dot(rotM,rotM.T),eye!(T)(ndim));
      if (err<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
      if (err2<T.mant_dig/4*3-11){
        Trace.formatln("err2: {}/{}",err2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
  
      auto err_3=minFeqrel2(rotM.T,rotM_2);
      if (err_3<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err_3,T.mant_dig);
        throw new Exception("inconsistent V1Ei EiV1",__FILE__,__LINE__);
      }
  }
  {
      auto rotM_2=rotateEiV!(NArray!(T,1),NArray!(T,2))(idx,n1,scale*eye!(T)(ndim));
      auto err_2=minFeqrel2(dot(rotM_2,n2),scale*n1);
      auto err2_2=minFeqrel2(dot(rotM_2,rotM_2.T),(scale*scale)*eye!(T)(ndim));
      if (err_2<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err_2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
      if (err2_2<T.mant_dig/4*3-11){
        Trace.formatln("err2: {}/{}",err2_2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
  
      auto rotM=rotateVEi!(NArray!(T,1),NArray!(T,2))(n1,idx,scale*eye!(T)(ndim));
      auto err=minFeqrel2(dot(rotM,n1),scale*n2);
      auto err2=minFeqrel2(dot(rotM,rotM.T),(scale*scale)*eye!(T)(ndim));
      if (err<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
      if (err2<T.mant_dig/4*3-11){
        Trace.formatln("err2: {}/{}",err2,T.mant_dig);
        throw new Exception("incorrect rotation",__FILE__,__LINE__);
      }
  }
}

/// invariant v1v2
void invariantRotVEi(T,int ndim,U=SizedRandomNArray!(T,ndim,ndim))(SizedRandomNArray!(T,ndim) v1,size_t idim,U m){
    size_t idx=idim%ndim;
    auto n1=v1.arr/norm2(v1.arr);

    auto rotM=rotateVEi!(NArray!(T,1),NArray!(T,2))(n1,idx,m.arr.dup);
    auto rotM2=rotateEiV!(NArray!(T,1),NArray!(T,2))(idx,n1,rotM);
    auto err=minFeqrel2(rotM2,m.arr);
    if (err<T.mant_dig/4*3-11){
        Trace.formatln("err: {}/{}",err,T.mant_dig);
        throw new Exception("non invariant",__FILE__,__LINE__);
    }
}

/// invariant v1v2
void sameVEiVV(T,int ndim,U=SizedRandomNArray!(T,ndim,ndim))(SizedRandomNArray!(T,ndim) v1,size_t idim,U m){
    size_t idx=idim%ndim;
    auto n1=v1.arr/norm2(v1.arr);
    auto n2=zeros!(T)(ndim);
    n2[idx]=cast(T)1;

    {
        auto rotM=rotateVEi!(NArray!(T,1),NArray!(T,2))(n1,idx,m.arr.dup);
        auto rotM2=rotateVV!(NArray!(T,1),NArray!(T,2))(n1,n2,m.arr.dup);
        auto err=minFeqrel2(rotM2,rotM);
        if (err<T.mant_dig/4*3-11){
            Trace.formatln("err: {}/{}",err,T.mant_dig);
            throw new Exception("VEi VV diff",__FILE__,__LINE__);
        }
    }
    {
        auto rotM=rotateEiV!(NArray!(T,1),NArray!(T,2))(idx,n1,m.arr.dup);
        auto rotM2=rotateVV!(NArray!(T,1),NArray!(T,2))(n2,n1,m.arr.dup);
        auto err=minFeqrel2(rotM2,rotM);
        if (err<T.mant_dig/4*3-11){
            Trace.formatln("err: {}/{}",err,T.mant_dig);
            throw new Exception("EiV VV diff",__FILE__,__LINE__);
        }
    }
}

void addRotTstToCollection(T,int ndim)(TestCollection coll){
    autoInitTst.testNoFailF("rotV1V2!("~T.stringof~")",&doTestRotV1V2!(T,ndim),
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("rotV1J!("~T.stringof~")",&doTestRotV1J!(T,ndim),
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("rotV1V2!("~T.stringof~")Invar",&invariantRotV1V2!(T,ndim),
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("rotV1J!("~T.stringof~")Invar",&invariantRotVEi!(T,ndim),
        __LINE__,__FILE__,coll);
    autoInitTst.testNoFailF("rotVV_V1J!("~T.stringof~")",&sameVEiVV!(T,ndim),
        __LINE__,__FILE__,coll);
}

/// rotation tests
TestCollection rotateTests(TestCollection superColl){
    TestCollection coll=new TestCollection("rotV1V2",
        __LINE__,__FILE__,superColl);
    addRotTstToCollection!(float,3)(coll);
    addRotTstToCollection!(double,3)(coll);
    addRotTstToCollection!(real,3)(coll);
    return coll;
}

/// collection of all the tests on the util modules
TestCollection utilTests(TestCollection superColl){
    TestCollection coll=new TestCollection("util",
        __LINE__,__FILE__,superColl);
    rotateTests(coll);
    return coll;
}
