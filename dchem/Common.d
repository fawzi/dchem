/// common declarations for dchem
module dchem.Common;
import blip.narray.NArray;
public import blip.narray.NArray: index_type;
import xf.omg.core.LinearAlgebra;

alias double Real;
alias zeros!(Real) zerosR;
alias ones!(Real) onesR;
alias empty!(Real) emptyR;
alias index_type[] IndexArray;

// small vectors (structs)
alias Vector!(Real, 2)  vec2R;
alias Vector!(Real, 3)  vec3R;
alias Vector!(Real, 4)  vec4R;

alias Quaternion!(Real) Quat;

alias Matrix!(Real, 2, 2)   Mat2;
alias Matrix!(Real, 3, 3)   Mat3;
alias Matrix!(Real, 3, 4)   Mat34;
alias Matrix!(Real, 4, 4)   Mat4;

/// converts a vector to a NArray (copy)
NArray!(T,1) vToNA(T,int dim_)(Vector!(T,dim_) v){
    auto res=empty!(T)([dim_]);
    *cast(typeof(v.tuple)*)res.startPtrArray=v.tuple;
    return res;
}
/// converts a matrix to a NArray (copy)
NArray!(T,2) m2NA(T,int rows_, int cols_)(Matrix!(T,rows_,cols_) m){
    auto res=empty!(T)([rows_,cols_],true);
    *cast(typeof(m.tuple)*)res.startPtrArray=m.tuple;
    return res;
}
/// converts a Quaternion to a NArray (copy)
NArray!(T,1) q2NA(T)(Quaternion!(T) q){
    auto res=empty!(T)([4]);
    res.startPtrArray[0]=q.x;
    res.startPtrArray[1]=q.y;
    res.startPtrArray[2]=q.z;
    res.startPtrArray[3]=q.w;
    return res;
}
/// converts a vector to a NArray (shares memory)
NArray!(T,1) vToNA(T,int dim_)(Vector!(T,dim_) *v){
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
