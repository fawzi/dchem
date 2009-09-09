module dchem.sys.Cell;

/// conversion of a,b,c,alpha,beta,gamma to h
NArray!(Real,2) cellParam2h(Real[] params){
    if (params.length==3){
        auto res=zeros!(Real)([3,3]);
        res.diag()=params;
        return res;
    } else if (params.length==6){
        return cellParam2h(params[0],params[1],params[2],params[3],params[4],params[5]);
    } else {
        throw new Exception("invalid number of params",__FILE__,__LINE__);
    }
}

/// conversion of a,b,c,alpha,beta,gamma to h
NArray!(Real,2) cellParam2h(Real a,Real b,Real c,Real alpha,Real beta,Real gamma){
    auto deg2rad=pi/180.0;
    auto cosBC=cos(deg2rad*alpha);
    auto sinBC=sin(deg2rad*alpha);
    auto cosAC=cos(deg2rad*beta);
    auto cosAB=cos(deg2rad*gamma);
    auto sinAB=sin(deg2rad*gamma);

    auto Ax=a;
    auto Bx=b*cosAB;
    auto By=b*sinAB;
    
    Real Cx,Cy,Cz;
    
    // If sinAB is zero, then we can't determine C uniquely since it's defined
    // in terms of the angle between A and B.
    if (sinAB!=0){
        Cx=cosAC;
        Cy=(cosBC - cosAC * cosAB) / sinAB;
        Cz=sqrt(1.0 - Cx*Cx - Cy*Cy);
    } else {
        Cx=0.0;
        Cy=0.0;
        Cz=0.0;
    }
    
    return a2na([[Ax,Bx,c*Cx],
        [0.0,By,c*Cy],
        [0.0,0.0,c*Cz]]);
}

/// finds the cell parameters (same units as h, angles in degrees) for a given cell
Real[] h2CellParam(NArray!(Real,2)h,Real[] params=null){
    auto rad2deg=180.0/pi;
    params.length=6;
    auto ang=dot(h,h.T);
    params[0]=sqrt(ang[0,0]);
    params[1]=sqrt(ang[1,1]);
    params[2]=sqrt(ang[2,2]);
    params[3]=atan2(ang[1,2],params[1]*params[2])*rad2deg;
    params[4]=atan2(ang[0,1],params[0]*params[1])*rad2deg;
    params[5]=atan2(ang[0,1],params[0]*params[1])*rad2deg;
    return params;
}

void testCell2Param(Real a,Real b,Real c,Real alpha,Real beta,Real gamma,GaussReal d0,GaussReal d1, GaussReal d2){
    auto h=cellParam2h(a,b,c,alpha,beta,gamma);
    auto param=h2CellParam(h);
    error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        Sdout("error:")(error).newline;
        throw new Exception("Error too big",__FILE__,__LINE__);
    }
    dir=
    auto param=h2CellParam(h);
    error=abs(param[0]-a)+abs(param[1]-b)+abs(param[2]-c)+abs(param[3]-alpha)
        +abs(param[4]-beta)+abs(param[5]-gamma);
    if (error>1.e-9) {
        Sdout("error2:")(error).newline;
        throw new Exception("Error2 too big",__FILE__,__LINE__);
    }
}

/++
 +  description of the simulation cell
 +   - periodicity: an array with 1 if that dimension is periodic, 0 if it isn't
 +   - h: matrix that maps the cell [0,1]^3 to the real cell in atomic units
 +   - h_inv: inverse of h
 +   - x0: shift of the origin (the mapping reduced points -> real points is
 +     r_real=dot(h,r_red)+x0
 +/
class Cell
{
    int[3] periodic;
    NArray!(Real,2) h,hInv;
    NArray!(Real,1) x0;

    this(NArray!(Real,2)h,int[3] periodic,NArray!(Real,2)hInv=null){
        this.h=h;
        this.periodic=periodic;
        if (hInv is null){
            this.hInv=LinAlg.inv(h);
        } else {
            this.hInv=hInv;
        }
    }
    Cell dup(){
        return new Cell(h.dup(),periodic,hInv.dup);
    }
}
