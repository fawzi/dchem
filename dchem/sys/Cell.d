/// cell of a system
module dchem.sys.Cell;
import dchem.Common;
import blip.serialization.Serialization;
import blip.narray.NArray;
import blip.serialization.StringSerialize;
import tango.math.Math;

/// conversion of a,b,c,alpha,beta,gamma to h
Matrix!(T, 3, 3) cellParam2h(T)(T a,T b,T c,T alpha,T beta,T gamma){
    auto deg2rad=PI/180.0;
    auto cosBC=cos(deg2rad*alpha);
    auto sinBC=sin(deg2rad*alpha);
    auto cosAC=cos(deg2rad*beta);
    auto cosAB=cos(deg2rad*gamma);
    auto sinAB=sin(deg2rad*gamma);

    auto Ax=a;
    auto Bx=b*cosAB;
    auto By=b*sinAB;
    
    T Cx,Cy,Cz;
    
    // If sinAB is zero, then we can't determine C uniquely since it's defined
    // in terms of the angle between A and B.
    if (sinAB!=0){
        Cx=cosAC;
        Cy=(cosBC - cosAC * cosAB) / sinAB;
        auto d=1.0 - Cx*Cx - Cy*Cy;
        if (d>0){
            Cz=sqrt(d);
        } else {
            Cz=0.0;
        }
    } else {
        Cx=0.0;
        Cy=0.0;
        Cz=0.0;
    }
    
    return Matrix!(T,3,3)([Ax,0,0,Bx,By,0,c*Cx,c*Cy,c*Cz]);
}
/// ditto
Matrix!(T,3,3) cellParamNArr2h(T)(NArray!(T,1) params){
    assert(params.shape[0]==6,"you should give 6 parameters");
    return cellParam2h(params[0],params[1],params[2],params[3],params[4],params[5]);
}
/// ditto
Matrix!(T,3,3) cellParamArr2h(T)(T[] params){
    assert(params.length==6,"you should give 6 parameters");
    return cellParam2h(params[0],params[1],params[2],params[3],params[4],params[5]);
}

/// finds the cell parameters (same units as h, angles in degrees) for a given cell
Real[] h2CellParam(NArray!(Real,2)h,Real[] params=null){
    auto rad2deg=180.0/PI;
    params.length=6;
    auto ang=dot(h.T,h);
    params[0]=sqrt(ang[0,0]);
    params[1]=sqrt(ang[1,1]);
    params[2]=sqrt(ang[2,2]);
    params[3]=acos(ang[1,2]/(params[1]*params[2]))*rad2deg;
    params[4]=acos(ang[0,2]/(params[0]*params[2]))*rad2deg;
    params[5]=acos(ang[0,1]/(params[0]*params[1]))*rad2deg;
    return params;
}

/// finds the cell parameters (same units as h, angles in degrees) for a given cell
Real[] h2CellParam(Matrix!(Real,3,3)m,Real[] params=null){
    return h2CellParam(m2NAC(m),params);
}


/// implicitly assumes changes only in a
char[] cellLoopMixin(char[][] names,char[]op){
    char[] res=`
    {
        void doVOp(`;
    foreach (i,n;names){
        if (i!=0) res~=",";
        res~="typeof("~n~".x0)"~n;
    }
    res~=`){`;
    res~=op;
    res~=`
        }
        doVOp(`;
    foreach (i,n;names){
        if (i!=0) res~=",";
        res~=n~".x0";
    }
    res~=`);
        void doMOp(`;
    foreach (i,n;names){
        if (i!=0) res~=",";
        res~="typeof("~n~".h)"~n;
    }
    res~=`){`;
    res~=op;
    res~=`
        }
        doMOp(`;
    foreach (i,n;names){
        if (i!=0) res~=",";
        res~=n~".h";
    }
    res~=`);
    }`;
    return res;
}

/++
 +  description of the simulation cell (this should normally behave as an immutable object)
 +   - periodicity: an array with 1 if that dimension is periodic, 0 if it isn't
 +   - h: matrix that maps the cell [0,1]^3 to the real cell in atomic units
 +   - h_inv: inverse of h
 +   - x0: shift of the origin (the mapping reduced points -> real points is
 +     r_real=dot(h,r_red)+x0
 +/
class Cell(T)
{
    alias T dtype;
    Matrix!(T,3,3) h,_hInv;
    Vector!(T,3) x0;
    int[3] periodic;
    bool hInvOk;

    this(){} // just for serialization
    
    this(Matrix!(T,3,3) h,int[3] periodic){
        this(h,periodic,Vector!(T,3).zero);
    }
    this(Matrix!(T,3,3) h,int[3] periodic,Vector!(T,3) x0,bool hInvOk=false){
        this.h=h;
        this.periodic[]=periodic;
        this._hInv=hInv;
        this.hInvOk=hInvOk;
        this.x0=x0;
    }
    
    Matrix!(T,3,3) hInv(){
        if (!hInvOk){
            synchronized(this){
                if (!hInvOk){
                    /+if (h==Matrix!(T,3,3).zero){
                        _hInv=Matrix!(T,3,3).identity;
                    } else {+/
                        _hInv=h.inverse;
                    /+}+/
                    hInvOk=true;
                }
            }
        }
        return _hInv;
    }
    
    static void bypaxOp(V)(ref Cell!(T) y, Cell!(V)x,T a=cscalar!(T,1),T b=cscalar!(T,1)){
        auto res=y.dup();
        res.opBypax(x,a,b);
        y=res;
    }
    Cell!(T) opMul(V)(V x){
        auto res=this.dup();
        static if (is(typeof(res.opMulAssign(x)))){
            res.opMulAssign(cast(T)x);
        } else {
            static assert(0,"cannot perform Cell!("~T.stringof~".opMul with "~V.stringof~")");
        }
        return res;
    }
    /// y.opBypax(x,a,b): y = ax+by
    private void opBypax(V)(Cell!(V) x,T a=cscalar!(T,1),T b=cscalar!(T,1)){
        auto y=this;
        h.opBypax(h,a,b);
        mixin(cellLoopMixin(["y","x"],"y.opBypax(x,a,b);"));
    }
    
    private void opMulAssign()(T scale){
        auto x=this;
        mixin(cellLoopMixin(["x"],"x*=scale;"));
    }
    private void opMulAssign(V)(Cell!(V) y){
        auto x=this;
        mixin(cellLoopMixin(["x","y"],"x *= y;"));
    }
    
    private void opSliceAssign(V)(ref Cell!(V) c2){
        periodic[]=c2.periodic;
        h.set(c2.h);
        _hInv.set(c2._hInv);// recalculate inverse in case the precision of T> precision of V?
        hInvOk=c2.hInvOk;
        x0.set(c2.x0);
    }
    private void opSliceAssign()(T val){
        h[]=val;
        hInvOk=false;
        x0[]=val;
    }
    Cell!(V) dupT(V=T)(){
        Cell!(V) res=new Cell!(V)();
        res[]=this;
        return res;
    }
    alias dupT!(T) dup;
    mixin(serializeSome("","periodic|h|x0"));
    mixin printOut!();
}

