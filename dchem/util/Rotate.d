/// rotations
module dchem.util.Rotate;
import tango.math.IEEE;
import blip.narray.NArray;

/// rotation from one *unit* vector v1 to a *unit* vector v2
M rotateVV(V,M)(V v1,V v2,M m)
in{
    alias typeof(dot(v1,v2)) S;
    auto err=feqrel2(norm2(v1),cast(S)1);
    assert(err>S.mant_dig/4*3-6,"v1 has to be a unit vector");
    auto err2=feqrel2(norm2(v2),cast(S)1);
    assert(err2>S.mant_dig/4*3-6,"v2 has to be a unit vector");
}
body{
    scope c=dot(v1,v2);
    alias typeof(c) S;
    scope v1m=dot(v1,m);
    scope v2Ortho=v2.dup;
    v2Ortho-=c*v1;
    scope v2OrthoM=dot(v2Ortho,m);
    outer(v2-v1,v1m,m,cast(S)1,cast(S)1);
    if (c<0){ // avoid bad numerics for almost anti collinear vectors
        S vOrtC=-1/(-c+1);
        outer(vOrtC*v2Ortho+v1,v2OrthoM,m,cast(S)1,cast(S)1);
    } else {
        S vOrtC=-1/(c+1);
        outer(vOrtC*v2Ortho-v1,v2OrthoM,m,cast(S)1,cast(S)1);
    }
    return m;
}

/// rotation from one delta i vector to a *unit* vector v2
M rotateEiV(V,M)(size_t i,V v2,M m)
in{
    alias typeof(dot(v2,v2)) S;
    auto err=feqrel2(norm2(v2),cast(S)1);
    assert(err>S.mant_dig/4*3-6,"v2 has to be a unit vector");
}
body{
    scope c=v2[i];
    alias typeof(c) S;
    scope v1m=m[i];
    static if(is(typeof(v1m.dup()))){
      v1m=v1m.dup();
    }
    scope v2Ortho=v2.dup;
    auto v2Val=v2Ortho[i];
    v2Ortho[i]=cast(S)0;
    scope v2OrthoM=dot(v2Ortho,m);
    v2Ortho[i]=v2Val-cast(S)1;
    outer(v2Ortho,v1m,m,cast(S)1,cast(S)1);
    if (c<0){
        S vOrtC=-1/(-c+cast(S)1);
        v2Ortho*=vOrtC;
        v2Ortho[i]=cast(S)1;
        outer(v2Ortho,v2OrthoM,m,cast(S)1,cast(S)1);
    } else {
        S vOrtC=-1/(c+cast(S)1);
        v2Ortho*=vOrtC;
        v2Ortho[i]=cast(S)(-1);
        outer(v2Ortho,v2OrthoM,m,cast(S)1,cast(S)1);
    }
    return m;
}

/// rotation from one *unit* vector v2 to the vector delta_i (transpose of the previous one)
M rotateVEi(V,M)(V v2,size_t i,M m)
in{
    alias typeof(dot(v2,v2)) S;
    auto err=feqrel2(norm2(v2),cast(S)1);
    assert(err>S.mant_dig/4*3-6,"v2 has to be a unit vector");
}
body{
    scope c=v2[i];
    alias typeof(c) S;
    
    scope v1m=m[i];
    static if(is(typeof(v1m.dup()))){
      v1m=v1m.dup();
    }
    scope v2Ortho=v2.dup;
    auto v2Val=v2[i];
    v2Ortho[i]=v2Ortho[i]-cast(S)1;
    scope v2M=dot(v2Ortho,m);
    static if(is(typeof(m[i]+=v2M))){
        m[i]+=v2M;
    } else {
        m[i]=m[i]+v2M;
    }
    if (c<0){
        S vOrtC=-1/(-c+cast(S)1);
        v2M*=vOrtC;
        v2Ortho[i]=cast(S)0;
        v2M.axpby(v1m,vOrtC*(-v2Val+1)+1);
        outer(v2Ortho,v2M,m,cast(S)1,cast(S)1);
    } else {
        S vOrtC=-1/(c+cast(S)1);
        v2M*=vOrtC;
        v2Ortho[i]=cast(S)0;
        v2M.axpby(v1m,vOrtC*(-v2Val+1)-1);
        outer(v2Ortho,v2M,m,cast(S)1,cast(S)1);
    }
    return m;
}


