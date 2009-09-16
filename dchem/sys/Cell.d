/// cell of a system (testing part should probably be moved away)
module dchem.sys.Cell;
import dchem.Common;
import blip.serialization.Serialization;
import blip.narray.NArray;
import blip.serialization.StringSerialize;
import tango.math.Math;
import tango.util.log.Trace;

/// conversion of a,b,c,alpha,beta,gamma to h
NArray!(Real,2) cellParam2h(Real a,Real b,Real c,Real alpha,Real beta,Real gamma){
    auto deg2rad=PI/180.0;
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
        auto d=1.0 - Cx*Cx - Cy*Cy;
        if (d>0){
            Cz=sqrt(d);
        } else {
            Trace.formatln("Cx {} Cy {}",Cx,Cy);
            Cz=0.0;
        }
    } else {
        Cx=0.0;
        Cy=0.0;
        Cz=0.0;
    }
    
    return a2NA([[Ax,Bx,c*Cx],
        [0.0,By,c*Cy],
        [0.0,0.0,c*Cz]]);
}

/// conversion of a,b,c,alpha,beta,gamma to h
NArray!(Real,2) cellParam2h(NArray!(Real,1) params){
    assert(params.shape[0]==6,"you should give 6 parameters");
    return cellParam2h(params[0],params[1],params[2],params[3],params[4],params[5]);
}

/// conversion of a,b,c,alpha,beta,gamma to h
NArray!(Real,2) cellParam2h(Real[] params){
    assert(params.length==6,"you should give 6 parameters");
    return cellParam2h(params[0],params[1],params[2],params[3],params[4],params[5]);
}


/// finds the cell parameters (same units as h, angles in degrees) for a given cell
Real[] h2CellParam(NArray!(Real,2)h,Real[] params=null){
    auto rad2deg=180.0/PI;
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

    this(){} // just for serialization
    
    this(NArray!(Real,2)h,int[3] periodic,NArray!(Real,1) x0=null,NArray!(Real,2)hInv=null){
        this.h=h.dup;
        this.periodic[]=periodic;
        if (hInv is null){
            this.hInv=inv(h);
        } else {
            this.hInv=hInv.dup;
        }
        this.x0=x0.dup;
        if (this.x0 is null){ // leave null?
          this.x0=zerosR(3);
        }
    }
    Cell dup(){
        return new Cell(h,periodic,x0,hInv);
    }
    mixin(serializeSome("","periodic|h|hInv|x0"));
    char[] toString(){
        return serializeToString(this);
        /+char[] res;
        mixin("res=serializeToString(this);");
        return res;+/
    }
}



