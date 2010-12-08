/// common declarations for dchem
module dchem.Common;
import blip.narray.NArray;
public import blip.narray.NArray: index_type;
import blip.omg.core.LinearAlgebra;
public import blip.omg.core.Algebra: scalar,cscalar;
public import blip.omg.core.LinearAlgebra:Vector,Quaternion,Matrix;
public import blip.Comp;

alias real   HighP;
alias double Real;
alias float  LowP;

alias zeros!(Real) zerosR;
alias ones!(Real) onesR;
alias empty!(Real) emptyR;
alias index_type[] IndexArray;

/// converts a vector to a NArray (copy)
NArray!(T,1) v2NAC(T,int dim_)(Vector!(T,dim_) v){
    auto res=empty!(T)([dim_]);
    (cast(T*)res.startPtrArray)[0..dim_]=v.cell;
    return res;
}
/// converts a matrix to a NArray (copy)
NArray!(T,2) m2NAC(T,int rows_, int cols_)(Matrix!(T,rows_,cols_) m){
    auto res=empty!(T)([rows_,cols_],true);
    (cast(T*)res.startPtrArray)[0..rows_*cols_]=m.cell;
    //*cast(typeof(m.tuple)*)res.startPtrArray=m.tuple;
    return res;
}
/// converts a Quaternion to a NArray (copy)
NArray!(T,1) q2NAC(T)(Quaternion!(T) q){
    auto res=empty!(T)([4]);
    res.startPtrArray[0]=q.x;
    res.startPtrArray[1]=q.y;
    res.startPtrArray[2]=q.z;
    res.startPtrArray[3]=q.w;
    return res;
}
/// converts a vector to a NArray (shares memory)
NArray!(T,1) v2NA(T,int dim_)(Vector!(T,dim_) *v){
    auto res=NArray!(T,1)([cast(index_type)T.sizeof], [cast(index_type)dim_], 0,
        v.cell, 0);
    return res;
}
/// converts a matrix to a NArray (shares memory)
NArray!(T,2) m2NA(T,int rows_, int cols_)(Matrix!(T,rows_,cols_) *m){
    auto res=NArray!(T,2)([cast(index_type)T.sizeof,cast(index_type)(rows_*T.sizeof)],
        [rows_,cols_], 0,m.cell, 0);
    return res;
}
/// converts a Quaternion to a NArray (shares memory)
NArray!(T,1) q2NA(T)(Quaternion!(T) *q){
    auto res=NArray!(T,1)([cast(index_type)T.sizeof], [cast(index_type)dim_], 0,
        q.xyzw.cell, 0);
    return res;
}

Matrix!(T,3,3) na2m33(T)(NArray!(T,2) m){
    assert(m.shape[0]==3 && m.shape[1]==3,"invalid matrix shape to convert to 3x3 matrix");
    return Matrix!(T,3,3)(m.data);
}

Vector!(T,3) na2v3(T)(NArray!(T,1) v){
    assert(v.shape[0]==3,"invalid vector shape to convert to Vector 3");
    return Vector!(T,3)(v[0],v[1],v[2]);
}
