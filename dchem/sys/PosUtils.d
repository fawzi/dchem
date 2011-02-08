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

/// square of the euclidean norm of the distance between c1 and c2
Real norm22Threshold(T)(T[]c1,T[]c2,Real threshold=Real.max){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    Real rAtt=0;
    if (threshold>0){
        for(size_t i=0;i<c1.length;++i){
            rAtt+=pow2(c2[i]-c1[i]);
        }
    }
    return rAtt;
}

// distances between the atoms of two configurations for segmented arrays
void configDists(W,U,V,T)(Cell!(W) cell,SegmentedArray!(U)c1,SegmentedArray!(V)c2,SegmentedArray!(T)dist){
    static assert(is(W==T),"incompatible cell");
    static assert(is(U==Vector!(T,3)),"c1 should contain vectors");
    static assert(is(V==Vector!(T,3)),"c2 should contain vectors");
    assert(c1.kRange==c2.kRange,"different kRange of configurations");
    assert(c1.kRange==dist.kRange,"different kRange of positions");
    CellFlag flag=cell.flags;
    foreach(k;dist.kRange.pLoop){
        auto c1Arr=c1[k];
        auto c2Arr=c2[k];
        auto distArr=dist[k];
        switch(flag){
        case CellFlag.NoCell:
            nConfigDistsKind(c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        case CellFlag.Ortho:
            oConfigDistsKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        case CellFlag.General:
            gConfigDistsKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData(),distArr.basicData());
            break;
        default:
            assert(0);
        }
    }
}
// distances between the atoms of a kind of two configurations
void configDistsKind(W,U,V,T)(Cell!(W) cell,BulkArray!(U)c1Arr,BulkArray!(V)c2Arr,BulkArray!(T)distArr){
    assert(c1Arr.length==c2Arr.length,"configurations need to have the same length");
    assert(distArr.length==c2Arr.length,"positions have unexpected length");
    switch(cell.flags){
    case CellFlag.NoCell:
        nConfigDistsKind(c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    case CellFlag.Ortho:
        oConfigDistsKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    case CellFlag.General:
        gConfigDistsKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData(),distArr.basicData);
        break;
    default:
        assert(0);
    }
}

/// optimized loop for the pp periodicity
private char[] oConfigDsLoopMixin(int pp){
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

private char[] allOConfigDsLoopMixin(){
    char[] res=`
    switch (pp){`;
    for (int i=0;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=oConfigDsLoopMixin(i);
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
void nConfigDistsKind(T)(T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    mixin(oConfigDsLoopMixin(0));
}

// distances between the atoms of two configurations in the orthorombic periodic case
void oConfigDistsKind(T)(T[3]h_diag,int periodic,T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic;
    mixin(allOConfigDsLoopMixin());
}

/// optimized loop for the general case
private char[] gConfigDsLoopMixin(int pp){
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

private char[] allGConfigDsLoopMixin(){
    char[] res=`
    switch (pp){
    case 0:`;
    res~=oConfigDsLoopMixin(0);
    res~=`
        break;`;
    for (int i=1;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=gConfigDsLoopMixin(i);
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
void gConfigDistsKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,T[]dist){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    auto pp=periodic;
    mixin(allGConfigDsLoopMixin());
}

// distances between the atoms of two configurations, automatically chooses the optimal solution
void configDistsKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,T[]dist,CellFlag flags=CellFlag.Unknown){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    assert(c1.length/3==dist.length,"the distances should be one third of the positions");
    if (flags==CellFlag.Unknown) flags=flagForCell(h,periodic);
    switch(flags){
    case CellFlag.NoCell:
        nConfigDistsKind(c1,c2,dist);
        break;
    case CellFlag.Ortho:
        oConfigDistsKind([h[0],h[4],h[8]],periodic,c1,c2,dist);
        break;
    case CellFlag.General:
        gConfigDistsKind(h,hInv,periodic,c1,c2,dist);
        break;
    default:
        assert(0);
    }
}

//----------------

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
    auto pp=periodic;
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
    auto pp=periodic;
    mixin(allGOneToNLoopMixin());
}

// distances between one atom and a set of others automatically choosing the optimal method
void oneToNDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,out T[]dist,CellFlag flags=CellFlag.Unknown){
    assert(c1.length==3,"c1 should have 3 components");
    assert(c2.length%3==0,"c2 should have a length that is a multiple of 3");
    assert(c2.length/3==dist.length,"the distances should be one third of the positions");
    if (flags==CellFlag.Unknown) flags=flagForCell(h,periodic);
    switch(flags){
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

//----------------

/// distances between the atoms of two configurations for segmented arrays
Real configDist(W,U,V,T=W)(Cell!(W) cell,SegmentedArray!(U)c1,SegmentedArray!(V)c2,Real threshold=Real.max){
    static assert(is(U==Vector!(W,3)),"c1 should contain vectors");
    static assert(is(V==Vector!(W,3)),"c2 should contain vectors");
    assert(c1.kRange==c2.kRange,"different kRange of configurations");
    CellFlag flag=cell.flags;
    Real res=0;
    foreach(k;c1.kRange.pLoop){
        auto c1Arr=c1[k];
        auto c2Arr=c2[k];
        Real dAtt;
        switch(flag){
        case CellFlag.NoCell:
            dAtt=nConfigDistKind(c1Arr.basicData(),c2Arr.basicData(),threshold-res);
            break;
        case CellFlag.Ortho:
            dAtt=oConfigDistKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData(),threshold-res);
            break;
        case CellFlag.General:
            dAtt=gConfigDistKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData(),threshold-res);
            break;
        default:
            assert(0);
        }
        atomicAdd(res,dAtt);
        if (res>threshold) break;
    }
    return res;
}
/// distances between the atoms of a kind of two configurations
Real configDistKind(W,U,V,T=W)(Cell!(W) cell,BulkArray!(U)c1Arr,BulkArray!(V)c2Arr,Real threshold=Real.max){
    assert(c1Arr.length==c2Arr.length,"configurations need to have the same length");
    Real res=0;
    switch(cell.flags){
    case CellFlag.NoCell:
        res=nConfigDistKind(c1Arr.basicData(),c2Arr.basicData(),threshold);
        break;
    case CellFlag.Ortho:
        res=oConfigDistKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData(),threshold);
        break;
    case CellFlag.General:
        res=gConfigDistKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData(),threshold);
        break;
    default:
        assert(0);
    }
    return res;
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
        auto lenAtt=c2.length/3;
        for(size_t i=0;i<lenAtt;++i){
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
            rAtt+=pow2(r1)+pow2(r2)+pow2(r3);
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
Real nConfigDistKind(T)(T[]c1,T[]c2,Real threshold){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    Real rAtt=0;
    mixin(oConfigDLoopMixin(0));
    return rAtt;
}

// distances between the atoms of two configurations in the orthorombic periodic case
Real oConfigDistKind(T)(T[3]h_diag,int periodic,T[]c1,T[]c2,T threshold){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    Real rAtt=0;
    mixin(allOConfigDLoopMixin());
    return rAtt;
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
        auto lenAtt=c2.length/3;
        for(size_t i=0;i<lenAtt;++i){
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
            rAtt+=pow2(p1)+pow2(p2)+pow2(p3);
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

/// distances between the atoms of two configurations in the general case (first image convention)
/// not that in some peculiar cases this is different from the closest image
Real gConfigDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,Real threshold){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    Real rAtt=0;
    mixin(allGConfigDLoopMixin());
    return rAtt;
}

/// distances between the atoms of two configurations, automatically chooses the optimal solution
Real configDistKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,Real threshold,CellFlag flags=CellFlag.Unknown){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    if (flags==CellFlag.Unknown) flags=flagForCell(h,periodic);
    switch(flags){
    case CellFlag.NoCell:
        return nConfigDistKind(c1,c2,dist,threshold);
        break;
    case CellFlag.Ortho:
        return oConfigDistKind([h[0],h[4],h[8]],periodic,c1,c2,threshold);
        break;
    case CellFlag.General:
        return gConfigDistKind(h,hInv,periodic,c1,c2,threshold);
        break;
    default:
        assert(0);
    }
}

//----------------

/// first image convention between the atoms of two configurations for segmented arrays
void makeClose(W,U,V)(Cell!(W) cell,SegmentedArray!(U)c1,SegmentedArray!(V)c2){
    static assert(is(W==T),"incompatible cell");
    static assert(is(U==Vector!(T,3)),"c1 should contain vectors");
    static assert(is(V==Vector!(T,3)),"c2 should contain vectors");
    assert(c1.kRange==c2.kRange,"different kRange of configurations");
    CellFlag flag=cell.flags;
    foreach(k;c1.kRange.pLoop){
        auto c1Arr=c1[k];
        auto c2Arr=c2[k];
        switch(flag){
        case CellFlag.NoCell:
            break;
        case CellFlag.Ortho:
            oMakeCloseKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData());
            break;
        case CellFlag.General:
            gMakeCloseKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
                c1Arr.basicData(),c2Arr.basicData());
            break;
        default:
            assert(0);
        }
    }
}
/// apply first image convention to c2 between the atoms of a kind of two configurations
void makeCloseKind(W,U,V)(Cell!(W) cell,BulkArray!(U)c1Arr,BulkArray!(V)c2Arr){
    assert(c1Arr.length==c2Arr.length,"configurations need to have the same length");
    switch(cell.flags){
    case CellFlag.NoCell:
        break;
    case CellFlag.Ortho:
        oMakeCloseKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData());
        break;
    case CellFlag.General:
        gMakeCloseKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,
            c1Arr.basicData(),c2Arr.basicData());
        break;
    default:
        assert(0);
    }
}

/// optimized loop for the pp periodicity
private char[] oMakeCloseLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        auto lenAtt=c2.length/3;
        for(size_t i=0;i<lenAtt;++i){
            auto i3=3*i;`;
    if ((pp&1)!=0){
        res~=`
            auto r1=(c2[i3]-c1[i3]);
            r1=r1/h_diag[0];
            r1=h_diag[0]*(r1-floor(r1+half));
            c2[i3+0]=c1[i3+0]+r1;`;
    }
    if ((pp&2)!=0){
        res~=`
            auto r2=(c2[i3+1]-c1[i3+1]);
            r2=r2/h_diag[1];
            r2=h_diag[1]*(r2-floor(r2+half));
            c2[i3+1]=c1[i3+1]+r2;`;
    }
    if ((pp&4)!=0){
        res~=`
            auto r3=(c2[i3+2]-c1[i3+2]);
            r3=r3/h_diag[2];
            r3=h_diag[2]*(r3-floor(r3+half));
            c2[i3+2]=c1[i3+2]+r3;`;
    }
    res~=`
        }
    }
    `;
    return res;
}

private char[] allOMakeCloseLoopMixin(){
    char[] res=`
    switch (pp){`;
    for (int i=0;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=oMakeCloseLoopMixin(i);
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

// applies the first image convention to c2 in the orthorombic periodic case
void oMakeCloseKind(T)(T[3]h_diag,int periodic,T[]c1,T[]c2){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    mixin(allOMakeCloseLoopMixin());
}

/// optimized loop for the general case
private char[] gMakeCloseLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        auto lenAtt=c2.length/3;
        for(size_t i=0;i<lenAtt;++i){
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
            c2[i3+0]=c1[i3+0]+h[0]*r1+h[1]*r2+h[2]*r3;
            c2[i3+1]=c1[i3+1]+h[3]*r1+h[4]*r2+h[5]*r3;
            c2[i3+2]=c1[i3+2]+h[6]*r1+h[7]*r2+h[8]*r3;
        }
    }
    `;
    return res;
}

private char[] allGMakeCloseLoopMixin(){
    char[] res=`
    switch (pp){
    case 0:
        break;`;
    for (int i=1;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=gMakeCloseLoopMixin(i);
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

/// applies first image convention to c2
void gMakeCloseKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    mixin(allGMakeCloseLoopMixin());
}

// makes c2 close to c1, automatically chooses the optimal solution
void makeCloseKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,T[]c2,CellFlag flags=CellFlag.Unknown){
    assert(c1.length==c2.length,"the two configurations should have equal length");
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    if (flags==CellFlag.Unknown) flags=flagForCell(h,periodic);
    switch(flags){
    case CellFlag.NoCell:
        break;
    case CellFlag.Ortho:
        oMakeCloseKind([h[0],h[4],h[8]],periodic,c1,c2);
        break;
    case CellFlag.General:
        gMakeCloseKind(h,hInv,periodic,c1,c2);
        break;
    default:
        assert(0);
    }
}

//----------------

/// wraps the vector deltaX so that it is inside the first cell
void wrap(W,U)(Cell!(W) cell,SegmentedArray!(U)deltaX){
    static assert(is(W==T),"incompatible cell");
    static assert(is(U==Vector!(T,3)),"c1 should contain vectors");
    CellFlag flag=cell.flags;
    foreach(k;c1.kRange.pLoop){
        auto c1Arr=c1[k];
        auto c2Arr=c2[k];
        switch(flag){
        case CellFlag.NoCell:
            nWrapKind(c1Arr.basicData());
            break;
        case CellFlag.Ortho:
            oWrapKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
                c1Arr.basicData());
            break;
        case CellFlag.General:
            gWrapKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,c1Arr.basicData());
            break;
        default:
            assert(0);
        }
    }
}
/// distances between the atoms of a kind of two configurations
void wrapKind(W,U)(Cell!(W) cell,BulkArray!(U)c1Arr){
    switch(cell.flags){
    case CellFlag.NoCell:
        nWrapKind(c1Arr.basicData());
        break;
    case CellFlag.Ortho:
        oWrapKind([cell.h.cgetRC!(0,0)(),cell.h.cgetRC!(1,1)(),cell.h.cgetRC!(2,2)()],cell.periodicFlags,
            c1Arr.basicData());
        break;
    case CellFlag.General:
        gWrapKind(cell.h.cell,cell.hInv.cell,cell.periodicFlags,c1Arr.basicData());
        break;
    default:
        assert(0);
    }
}

/// optimized loop for the pp periodicity
private char[] oWrapLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        auto lenAtt=c1.length/3;
        for(size_t i=0;i<lenAtt;++i){
            auto i3=3*i;`;
    if ((pp&1)!=0){
        res~=`
            auto r1=c1[i3];
            r1=r1/h_diag[0];
            r1=h_diag[0]*(r1-floor(r1+half));
            c1[i3+0]=r1;`;
    }
    if ((pp&2)!=0){
        res~=`
            auto r2=c1[i3+1];
            r2=r2/h_diag[1];
            r2=h_diag[1]*(r2-floor(r2+half));
            c1[i3+1]=r2;`;
    }
    if ((pp&4)!=0){
        res~=`
            auto r3=c1[i3+2];
            r3=r3/h_diag[2];
            r3=h_diag[2]*(r3-floor(r3+half));
            c1[i3+2]=r3;`;
    }
    res~=`
        }
    }
    `;
    return res;
}

private char[] allOWrapLoopMixin(){
    char[] res=`
    switch (pp){`;
    for (int i=0;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=oWrapLoopMixin(i);
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

// distances between the atoms of two configurations in the orthorombic periodic case
void oWrapKind(T)(T[3]h_diag,int periodic,T[]c1){
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    mixin(allOWrapLoopMixin());
}

/// optimized loop for the general case
private char[] gWrapLoopMixin(int pp){
    char[] res=`
    {
        version(D_Version2){
            mixin("enum { half=(cast(T)1)/(cast(T)2) }");
        } else {
            const T half=(cast(T)1)/(cast(T)2);
        }
        auto lenAtt=c2.length/3;
        for(size_t i=0;i<lenAtt;++i){
            auto i3=3*i;
            auto p1=c1[i3];
            auto p2=c1[i3+1];
            auto p3=c1[i3+2];
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
            c1[i3+0]=h[0]*r1+h[1]*r2+h[2]*r3;
            c1[i3+1]=h[3]*r1+h[4]*r2+h[5]*r3;
            c1[i3+2]=h[6]*r1+h[7]*r2+h[8]*r3;
        }
    }
    `;
    return res;
}

private char[] allGWrapLoopMixin(){
    char[] res=`
    switch (pp){
    case 0:
        break;`;
    for (int i=1;i<8;++i){
        res~=`
    case `~ctfe_i2a(i)~`:`;
        res~=gWrapLoopMixin(i);
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
void gWrapKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1){
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    auto pp=periodic;
    mixin(allGWrapLoopMixin());
}

// distances between the atoms of two configurations, automatically chooses the optimal solution
void wrapKind(T)(T[9]h,T[9]hInv,int periodic,T[]c1,CellFlag flags=CellFlag.Unknown){
    assert(c1.length%3==0,"the configurations should have a length that is a multiple of 3");
    if (flags==CellFlag.Unknown) flags=flagForCell(h,periodic);
    switch(flags){
    case CellFlag.NoCell:
        break;
    case CellFlag.Ortho:
        oWrapKind([h[0],h[4],h[8]],periodic,c1);
        break;
    case CellFlag.General:
        gWrapKind(h,hInv,periodic,c1);
        break;
    default:
        assert(0);
    }
}

//----------------

version(InstantiateSome){
    private {
        // allocate something for testing purposes
        alias configDists!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyA;
        alias configDistsKind!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyB;
        alias configDistsKind!(Real) dummyC;
        alias configDist!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyD;
        alias configDistKind!(Real,Vector!(Real,3),Vector!(Real,3),Real) dummyE;
        alias configDistKind!(Real) dummyF;
        alias oneToNDistKind!(Real) dummyG;
    }
}