/// rotations
module dchem.util.Rotate;
import blip.narray.NArray;
import tango.util.log.Trace;
import tango.io.Stdout;

/// rotation from one *unit* vector v1 to a *unit* vector v2
M rotateVV(V,M)(V v1,V v2,M m){
    auto c=dot(v1,v2);
    alias typeof(c) S;
    auto s2=cast(S)1-c*c;
    S vOrtC;
    if (s2>0){
      vOrtC=(c-cast(S)1)/s2;
    } else {
      vOrtC=cast(S)(-1)/cast(S)2;
    }
    auto v1m=dot(v1,m);
    auto v2Ortho=v2.dup;
    v2Ortho-=c*v1;
    auto v2OrthoM=dot(v2Ortho,m);
    outer(v2-v1,v1m,m,cast(S)1,cast(S)1);
    outer(vOrtC*v2Ortho-v1,v2OrthoM,m,cast(S)1,cast(S)1);
    return m;
}

/// rotation from one delta i vector to a *unit* vector v2
M rotateEiV(V,M)(size_t i,V v2,M m){
    auto c=v2[i];
    alias typeof(c) S;
    auto s2=cast(S)1-c*c;
    S vOrtC;
    if (s2>0){
      vOrtC=(c-cast(S)1)/s2;
    } else {
      vOrtC=cast(S)(-1)/cast(S)2;
    }
    auto v1m=m[i];
    static if(is(typeof(v1m.dup()))){
      v1m=v1m.dup();
    }
    auto v2Ortho=v2.dup;
    auto v2OrthoM=dot(v2,m);
    v2OrthoM[i]=cast(S)0;
    v2Ortho[i]=v2Ortho[i]-cast(S)1;
    outer(v2Ortho,v1m,m,cast(S)1,cast(S)1);
    v2Ortho*=vOrtC;
    v2Ortho[i]=cast(S)(-1);
    outer(v2Ortho,v2OrthoM,m,cast(S)1,cast(S)1);
    return m;
}

/// rotation from one *unit* vector v2 to the vector delta_i (transpose of the previous one)
M rotateVEi(V,M)(V v2,size_t i,M m){
    auto c=v2[i];
    alias typeof(c) S;
    auto s2=cast(S)1-c*c;
    S vOrtC;
    if (s2>0){
      vOrtC=(c-cast(S)1)/s2;
    } else {
      vOrtC=cast(S)(-1)/cast(S)2;
    }
    auto v1m=m[i];
    static if(is(typeof(v1m.dup()))){
      v1m=v1m.dup();
    }
    auto v2Ortho=v2.dup;
    auto v2OrthoM=dot(v2,m);
    v2OrthoM[i]=cast(S)0;
    v2Ortho[i]=v2Ortho[i]-cast(S)1;
    outer(v1m,v2Ortho,m,cast(S)1,cast(S)1);
    v2Ortho*=vOrtC;
    v2Ortho[i]=cast(S)-1;
    outer(v2OrthoM,v2Ortho,m,cast(S)1,cast(S)1);
    return m;
}


