/// rotations
module dchem.Rotate;

/// rotation from one unit vector v1 to a unit vector v2
M rotate(M,V)(V v1,V v2,M m){
    real c,s;
    c=dot(v1,v2);
    res=outer(2*c*v2-v1,dot(v2,m))+outer(v1,dot(v2,m));
    if (c==1.0){
        s=1.0;
    } else {
        s=sqrt(1.0-c);
    }
    m-=outer(v1,dot(v1/s**2,m))-outer(v2-c*v1,dot(v2-c*v1,m)/s**2);
}

