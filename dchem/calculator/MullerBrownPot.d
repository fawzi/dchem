module dchem.calculator.MullerBrownPot;
import dchem.calculator.Calculator;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.sys.SegmentedArray;
import dchem.Common;
import dchem.input.RootInput;
import blip.io.Console;
import blip.container.GrowableArray;
import blip.io.BasicIO;
import dchem.calculator.ProcContext;
import blip.serialization.Serialization;
import dchem.input.ReadIn2PSys: artificialPSys;
import blip.t.math.Math;

const Real[4] MBPreFactor = [-200.0,-100.0,-170.0,15.0];
const Real[4] MBa = [-1.0,-1.0,-6.5,0.7];
const Real[4] MBb = [0.0,0.0,11.0,0.6];
const Real[4] MBc = [-10.0,-10.0,-6.5,0.7];
const Real[4] MBx0 = [1.0,0.0,-0.5,-1.0];
const Real[4] MBy0 = [0.0,0.5,1.5,1.0];

class MullerBrownPot: Method {
    Real[] a=MBa,b=MBb,c=MBc,x0=MBx0,y0=MBy0,preFactor=MBPreFactor;
    
    this(){}
    
    bool verify(CharSink sink){
        bool res=true;
        auto s=dumper(sink);
        if (a.length!=b.length || a.length!=c.length || a.length!=x0.length || a.length!=x0.length
            || a.length!=preFactor.length){
                s("Inconsistent lengths arrays of "~myFieldName);
            }
        return res;
    }
    
    void deriv(Real x,Real y,out Real dx, out Real dy){
        sout("deriv\n");
        auto len=a.length;
        dx=0;
        dy=0;
        sout("arrs\n");
        sout("a:")(a.length)(" b:")(b.length)(" c:")(c.length)(" x0:")(x0.length)(" y0:")(y0.length)(" preFactor:")(preFactor.length)("\n");
        for (size_t ifact=0;ifact<len;++ifact){
            dx=preFactor[ifact]*exp(a[ifact]*pow2((x - x0[ifact])) + b[ifact]*(x - x0[ifact])*(y - y0[ifact]) + c[ifact]*pow2(y - y0[ifact]))*(a[ifact]*2.0*(x - x0[ifact]) + b[ifact]*(y - y0[ifact]));
            dy=preFactor[ifact]*exp(a[ifact]*pow2((x - x0[ifact])) + b[ifact]*(x - x0[ifact])*(y - y0[ifact]) + c[ifact]*pow2(y - y0[ifact]))*(b[ifact]*(x - x0[ifact]) + c[ifact]*2.0*(y - y0[ifact]));
        }
        sout("derivEnd\n");
    }
    
    Real eval(Real x,Real y){
        Real res=0;
        auto len=a.length;
        for (size_t ifact=0;ifact<len;++ifact){
            res+=preFactor[ifact]*exp(a[ifact]*pow2((x - x0[ifact])) + b[ifact]*(x - x0[ifact])*(y - y0[ifact]) + c[ifact]*pow2(y - y0[ifact]));
        }
        return res;
    }
    /// returns a new context
    CalcContext getCalculator(bool wait,ubyte[]history){
        auto ctx=new MullerBrownPotContext(this,collectAppender(delegate void(CharSink s){
            dumper(s)("MBPotCtx")(ProcContext.instance.id)("-")(ProcContext.instance.localId.next());
        }));
        return ctx;
    }
    /// drops the history with the given id
    void dropHistory(ubyte[]history) { }
    /// clears all history
    void clearHistory() { }
    /// url to access this from other processes
    char[] exportedUrl(){
        assert(0,"to do");
    }
    
    // serialization stuff
    mixin(serializeSome("dchem.MullerBrownPot",
    `a|b|c|x0|y0|preFactor`));
    mixin printOut!();
    mixin myFieldMixin!();
}

class MullerBrownPotContext:CalcContext{
    MullerBrownPot pot;
    
    this(MullerBrownPot pot,char[] contextId){
        super(contextId);
        this.pot=pot;
        _pSysReal=artificialPSys!(Real)(0,0,2);
    }

    Method method(){
        return pot;
    }
    
    /// should collect the newly calculated energy
    void updateEF(bool updateE=true,bool updateF=true){
        auto sys=LocalPIndex(0,0);
        auto x=pSysReal.dynVars.x.dof[sys,0];
        auto y=pSysReal.dynVars.x.dof[sys,1];
        if (updateF){
            pSysReal.checkMddx();
            Real dx,dy;
            pot.deriv(x,y,dx,dy);
            pSysReal.dynVars.mddx.dof[sys,0]= -dx;
            pSysReal.dynVars.mddx.dof[sys,1]= -dy;
        }
        if (updateE){
            pSysReal.dynVars.potentialEnergy=pot.eval(x,y);
        }
    }
    
}