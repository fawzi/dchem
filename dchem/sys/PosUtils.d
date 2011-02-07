/// utilities for distances and positions of 3D vectors
///
/// could be easily wekened wrt. to acceptance of different types..., do it??
module dchem.sys.PosUtils;
import dchem.Common;
import blip.serialization.Serialization;
import blip.serialization.StringSerialize;
import blip.math.Math;
import dchem.sys.Cell;
import dchem.sys.SegmentedArray;
import blip.container.BulkArray;

// distances between the atoms of two configurations for segmented arrays
void configDist(W,U,V,T)(Cell!(W) cell,SegmentedArray!(U)c1,SegmentedArray!(V)c2,SegmentedArray!(T)dist){
    static assert(is(W==T),"incompatible cell");
    static assert(is(U==Vector!(T,3)),"c1 should contain vectors");
    static assert(is(V==Vector!(T,3)),"c2 should contain vectors");
    assert(c1.kRange==c2.kRange,"different kRange of configurations");
    assert(c1.kRange==dist.kRange,"different kRange of positions");
    CellFlag flag=flagForCell(cell.h.cell,cell.periodic);
    foreach(k;dist.kRange.pLoop){
        auto c1Arr=c1[k];
        auto c2Arr=c2[k];
        auto distArr=dist[k];
        switch(flag){
        case CellFlag.NoCell:
            nConfigDistKind(c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        case CellFlag.Ortho:
            oConfigDistKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodic,
                c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        case CellFlag.General:
            gConfigDistKind(cell.h.cell,cell.hInv.cell,cell.periodic,
                c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        default:
            assert(0);
        }
    }
}
// distances between the atoms of a kind of two configurations
void configDistKind(W,U,V,T)(Cell!(W) cell,BulkArray!(U)c1Arr,BulkArray!(V)c2Arr,BulkArray!(T)distArr){
    assert(c1Arr.length==c2Arr.length,"configurations need to have the same length");
    assert(distArr.length==c2Arr.length,"positions have unexpected length");
    switch(flagForCell(cell.h.cell,cell.periodic)){
    case CellFlag.NoCell:
        nConfigDistKind(c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    case CellFlag.Ortho:
        oConfigDistKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodic,
            c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    case CellFlag.General:
        gConfigDistKind(cell.h.cell,cell.hInv.cell,cell.periodic,
            c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    default:
        assert(0);
    }
}

/// optimized loop for the pp periodicity
private char[] oConfigDLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        for(size_t i=0;i<dist.length;++i){
            auto i3=3*i;
            auto r1=(c2[i3]-c1[i3]);`;
    if ((pp&1)!=0){
        res~=`
            r1=r1/h_diag[0];
            r1=h_diag[0]*(r1-floor(r1+half));`;
    }
    res~=`
            auto r2=(c2[i3+1]-c1[i3+1]);`;
    if ((pp&2)!=0){
        res~=`
            r2=r2/h_diag[1];
            r2=h_diag[1]*(r2-floor(r2+half));`;
    }
    res~=`
            auto r3=(c2[i3+2]-c1[i3+2]);`;
    if ((pp&4)!=0){
        res~=`
            r3=r3/h_diag[2];
            r3=h_diag[2]*(r3-floor(r3+half));`;
    }
    res~=`
            dist[i]=pow2(r1)+pow2(r2)+pow2(r3);
        }
    }
    `;
    return res;
}

private char[] allOConfigDLoopMixin(){
    char[] res=`
    switch (pp){`;
    for (int i=0;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=oConfigDLoopMixin(i);
        res~=`
        break;`;
    }
    res~=`
    default:
        assert(0);
    }
    `;
    return res;
}

// distances between the atoms of two configurations in the non periodic case
void nConfigDistKind(T)(T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    mixin(oConfigDLoopMixin(0));
}

// distances between the atoms of two configurations in the orthorombic periodic case
void oConfigDistKind(T)(T[3]h_diag,int periodic,T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic[0]+2*periodic[1]+4*periodic[2];
    mixin(allOConfigDLoopMixin());
}

/// optimized loop for the general case
private char[] gConfigDLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        for(size_t i=0;i<dist.length;++i){
            auto i3=3*i;
            auto p1=(c2[i3]-c1[i3]);
            auto p2=(c2[i3+1]-c1[i3+1]);
            auto p3=(c2[i3+2]-c1[i3+2]);
            auto r1=hInv[0]*p1+hInv[1]*p2+hInv[2]*p3;
            auto r2=hInv[3]*p1+hInv[4]*p2+hInv[5]*p3;
            auto r3=hInv[6]*p1+hInv[7]*p2+hInv[8]*p3;
            `;
    if ((pp&1)!=0){
        res~=`
            r1=(r1-floor(r1+half));`;
    }
    if ((pp&2)!=0){
        res~=`
            r2=(r2-floor(r2+half));`;
    }
    if ((pp&4)!=0){
        res~=`
            r3=(r3-floor(r3+half));`;
    }
    res~=`
            p1=h[0]*r1+h[1]*r2+h[2]*r3;
            p2=h[3]*r1+h[4]*r2+h[5]*r3;
            p3=h[6]*r1+h[7]*r2+h[8]*r3;
            dist[i]=pow2(p1)+pow2(p2)+pow2(p3);
        }
    }
    `;
    return res;
}

private char[] allGConfigDLoopMixin(){
    char[] res=`
    switch (pp){
    case 0:`;
    res~=oConfigDLoopMixin(0);
    res~=`
        break;`;
    for (int i=1;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=gConfigDLoopMixin(i);
        res~=`
        break;`;
    }
    res~=`
    default:
        assert(0);
    }
    `;
    return res;
}

// distances between the atoms of two configurations in the general case (first image convention)
// not that in some peculiar cases this is different from the closest image
void gConfigDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic[0]+2*periodic[1]+4*periodic[2];
    mixin(allGConfigDLoopMixin());
}

// distances between the atoms of two configurations, automatically chooses the optimal solution
void configDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    switch(flagForCell(h,periodic)){
    case CellFlag.NoCell:
        nConfigDistKind(c1,c2,dist);
        break;
    case CellFlag.Ortho:
        oConfigDistKind([h[0],h[4],h[8]],periodic,c1,c2,dist);
        break;
    case CellFlag.General:
        gConfigDistKind(h,hInv,periodic,c1,c2,dist);
        break;
    default:
        assert(0);
    }
}

/// optimized loop for the orthorombic case
private char[] oOneToNLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            enum { half=(cast(T)1)/(cast(T)2) }
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        for(size_t i=0;i<dist.length;++i){
            auto i3=3*i;
            auto r1=(c2[i3]-c1[0]);`;
    if ((pp&1)!=0){
        res~=`
            r1=r1/h_diag[0];
            r1=h_diag[0]*(r1-floor(r1+half));`;
    }
    res~=`
            auto r2=(c2[i3+1]-c1[1]);`;
    if ((pp&2)!=0){
        res~=`
            r2=r2/h_diag[1];
            r2=h_diag[1]*(r2-floor(r2+half));`;
    }
    res~=`
            auto r3=(c2[i3+2]-c1[2]);`;
    if ((pp&4)!=0){
        res~=`
            r3=r3/h_diag[2];
            r3=h_diag[2]*(r3-floor(r3+half));`;
    }
    res~=`
            dist[i]=pow2(r1)+pow2(r2)+pow2(r3);
        }
    }
    `;
    return res;
}

private char[] allOOneToNLoopMixin(){
    char[] res=`
    switch (pp){`;
    for (int i=0;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=oOneToNLoopMixin(i);
        res~=`
        break;`;
    }
    res~=`
    default:
        assert(0);
    }
    `;
    return res;
}

// distances between the atoms of two configurations in the non periodic case
void nOneToNDistKind(T)(T[]c1,T[]c2,T[]dist){
    assert(c1.length==3,"c1 should have 3 components");
    assert(c2.length%3==0,"c2 should have a length that is a multiple of 3");
    assert(c2.length/3==dist.length,"the distances should be one third of the positions");
    mixin(oOneToNLoopMixin(0));
}

// distances between the atoms of two configurations in the orthorombic periodic case
void oOneToNDistKind(T)(T[3]h_diag,int periodic,T[]c1,T[]c2,out T[]dist){
    assert(c1.length==3,"c1 should have 3 components");
    assert(c2.length%3==0,"c2 should have a length that is a multiple of 3");
    assert(c2.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic[0]+2*periodic[1]+4*periodic[2];
    mixin(allOOneToNLoopMixin());
}

/// optimized loop for 
private char[] gOneToNLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        for(size_t i=0;i<dist.length;++i){
            auto i3=3*i;
            auto p1=(c2[i3+0]-c1[0]);
            auto p2=(c2[i3+1]-c1[1]);
            auto p3=(c2[i3+2]-c1[2]);
            auto r1=hInv[0]*p1+hInv[1]*p2+hInv[2]*p3;
            auto r2=hInv[3]*p1+hInv[4]*p2+hInv[5]*p3;
            auto r3=hInv[6]*p1+hInv[7]*p2+hInv[8]*p3;
            `;
    if ((pp&1)!=0){
        res~=`
            r1=(r1-floor(r1+half));`;
    }
    if ((pp&2)!=0){
        res~=`
            r2=(r2-floor(r2+half));`;
    }
    if ((pp&4)!=0){
        res~=`
            r3=(r3-floor(r3+half));`;
    }
    res~=`
            p1=h[0]*r1+h[1]*r2+h[2]*r3;
            p2=h[3]*r1+h[4]*r2+h[5]*r3;
            p3=h[6]*r1+h[7]*r2+h[8]*r3;
            dist[i]=pow2(p1)+pow2(p2)+pow2(p3);
        }
    }
    `;
    return res;
}

private char[] allGOneToNLoopMixin(){
    char[] res=`
    switch (pp){
    case 0:`;
    res~=oOneToNLoopMixin(0);
    res~=`
        break;`;
    for (int i=1;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=gOneToNLoopMixin(i);
        res~=`
        break;`;
    }
    res~=`
    default:
        assert(0);
    }
    `;
    return res;
}

// distances between one atom and a set of others in the non orthorombic periodic case
void gOneToNDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,out T[]dist){
    assert(c1.length==3,"c1 should have 3 components");
    assert(c2.length%3==0,"c2 should have a length that is a multiple of 3");
    assert(c2.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic[0]+2*periodic[1]+4*periodic[2];
    mixin(allGOneToNLoopMixin());
}

// distances between one atom and a set of others automatically choosing the optimal method
void oneToNDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,out T[]dist){
    assert(c1.length==3,"c1 should have 3 components");
    assert(c2.length%3==0,"c2 should have a length that is a multiple of 3");
    assert(c2.length/3==dist.length,"the distances should be one third of the positions");
    switch(flagForCell(h,periodic)){
    case CellFlag.NoCell:
        nOneToNDistKind(c1,c2,dist);
        break;
    case CellFlag.Ortho:
        oOneToNDistKind([h[0],h[4],h[8]],periodic,c1,c2,dist);
        break;
    case CellFlag.General:
        gOneToNDistKind(h,hInv,periodic,c1,c2,dist);
        break;
    default:
        assert(0);
    }
}

version(InstantiateSome){
    private {
        // allocate something for testing purposes
        alias configDist!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyA;
        alias configDistKind!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyB;
        alias configDistKind!(Real) dummyC;
        alias oneToNDistKind!(Real) dummyD;
    }
}