/// constraint that removes translations and rotations by fixing 3 particles
module dchem.sys.RotTransConstraint;
import dchem.Common;
import blip.serialization.Serialization;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import blip.BasicModels;
import blip.container.BulkArray;
import blip.io.BasicIO;
import blip.io.Console:serr;
import blip.container.GrowableArray;
import tango.math.Math:max;
import dchem.input.RootInput;
import dchem.sys.DynVars;
import dchem.sys.ParticleRange;
import dchem.sys.Constraints;
import dchem.sys.SegmentedArray;
import blip.io.Console;
import blip.math.Math;
import dchem.util.Rotate;

struct RotTTransf(T){
    Matrix!(T,3,3) rot;
    Vector!(T,3) transl;
    
    Real apply(BulkArray!(Vector!(T,3))arr){
        Real err=0;
        foreach(ref v;arr){
            auto res=rot.opVecMul(v+transl);
            err+=(res-v).norm22;
            v=res;
        }
        return err;
    }
    Real apply(SegmentedArray!(Vector!(T,3))sarr){
        Real err=0;
        foreach(k;sarr.kRange){
            foreach(ref v;sarr[k]){
                auto res=rot.opVecMul(v+transl);
                err+=(res-v).norm22;
                v=res;
            }
        }
        return err;
    }
    Real addDiff(BulkArray!(Vector!(T,3))arr,BulkArray!(Vector!(T,3))arr2){
        Real err=0;
        foreach(i,v;arr){
            auto diff=rot.opVecMul(v+transl)-v;
            arr2[i]+=diff;
            err+=diff.norm22;
        }
        return err;
    }
    Real addDiff(SegmentedArray!(Vector!(T,3))sarr,
                   SegmentedArray!(Vector!(T,3))sarr2){
        Real err=0;
        foreach(k;sarr.kRange){
            auto arr2=sarr2[k];
            foreach(i,v;sarr[k]){
                auto diff=rot.opVecMul(v+transl)-v;
                err+=diff.norm22;
                arr2[i]+=diff;
            }
        }
        return err;
    }
    Real diff(BulkArray!(Vector!(T,3))arr){
        Real err=0;
        foreach(i,v;arr){
            auto diff=rot.opVecMul(v+transl)-v;
            err+=diff.norm22;
        }
        return err;
    }
    Real diff(SegmentedArray!(Vector!(T,3))sarr){
        Real err=0;
        foreach(k;sarr.kRange){
            foreach(i,v;sarr[k]){
                auto diff=rot.opVecMul(v+transl)-v;
                err+=diff.norm22;
            }
        }
        return err;
    }
}

class RotTransConstraintGen:ConstraintGen {
    InputField particles;

    mixin myFieldMixin!();
    mixin(serializeSome("dchem.RotTransConstraint",`
    particles: the three or more particles used to remove translations and rotations
    `));
    mixin printOut!();
    
    this(){}
    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        if (particles is null) {
            log("particles should be a ParticleRange containing the (at least) 3 particles used to remove rotations and translations in field ")(myFieldName)("\n");
            res=false;
        } else if ((cast(ParticleRange)(particles.contentObj())) is null){
             log("particles should be a ParticleRange containing the (at least) 3 particles used to remove rotations and translations and not ")
                 (particles.contentObj.classinfo.name)(" in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    ConstraintI!(Real) constraintReal(ParticleSys!(Real)templatePos){
        return constraintT!(Real)(templatePos);
    }
    ConstraintI!(LowP) constraintLowP(ParticleSys!(LowP)templatePos){
        return constraintT!(LowP)(templatePos);
    }
    
    RotTransConstraint!(T) constraintT(T)(ParticleSys!(T)templatePos){
        return new RotTransConstraint!(T)(this,templatePos);
    }
}

class RotTransConstraint(T): ConstraintI!(T){
    RotTransConstraintGen _constraintGen;
    
    this(RotTransConstraintGen input,ParticleSys!(T)templatePos){
        this._constraintGen=input;
    }
    
    RotTTransf!(T) rotTForPSys(ParticleSys!(T)pSys){
        RotTTransf!(T) res;
        res.rot=Matrix!(T,3,3).identity;
        res.transl=Vector!(T,3).zero;
        size_t i=0;
        auto pos=pSys.dynVars.x.pos;
        bool stop=false;
        foreach (pks;_constraintGen.particles.contentT!(ParticleRange)().loopOn(pSys.sysStruct)){
            if (stop) break;
            foreach(p;pks){
                if (stop) break;
                switch(i){
                case 0:
                    res.transl= -pos[p,0];
                    break;
                case 1:
                    auto n=(pos[p,0]+res.transl).normalized();
                    assert(abs(n.norm2-1)<1.e-5);
                    res.rot=rotateVEi(n,0,res.rot);
                    break;
                default:
                    auto nv=res.rot*pos[p,0];
                    nv.x=0;
                    if (nv.norm22>0){
                        nv.normalize();
                        res.rot=rotateVEi(nv,1,res.rot);
                        stop=true;
                    }
                    break;
                }
                ++i;
            }
        }
        if (stop) return res;
        if (pos[pSys.sysStruct.levels[0]].length>3){
            if (i<3){
                throw new Exception("RotTransConstraint of field "~_constraintGen.myFieldName~
                                   " has less than 3 particles to remove translations and rotations",
                                   __FILE__,__LINE__);
            } else {
                // perfect alignment
                sout("Warning: perfect alignment for RotTransConstraint of field "~_constraintGen.myFieldName~"\n");
            }
        }
        return res;
    }
    
    ConstraintGen constraintGen(){
        return _constraintGen;
    }

    Real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null){
        auto op=rotTForPSys(state);
        if (deriv!is null){
            return op.addDiff(state.dynVars.x.pos,deriv.dynVars.dx.pos);
        } else {
            return op.diff(state.dynVars.x.pos);
        }
    }
    
    Real applyR(ParticleSys!(T) state,bool partial=false){
        auto op=rotTForPSys(state);
        return op.apply(state.dynVars.x.pos);
    }
    
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null){
        Real a=0,b=0,c=0,d=0;
        long n=0;
        foreach(i,pk,lk,v;state.dynVars.dx.pos.sLoop()){
            a+=v.x-v.y;
            b+=v.x-v.z;
            c+=v.y-v.z;
            d+=v.x+v.y+v.z;
            ++n;
        }
        a/=2*n;
        b/=2*n;
        c/=2*n;
        d/=3*n;
        foreach(i,pk,lk,ref v;state.dynVars.dx.pos.sLoop()){
            v.x+=-a-b-d;
            v.y+=a-c-d;
            v.z+=b+c-d;
        }
        if (diff!is null){
            foreach(i,pk,lk,ref v;diff.dynVars.dx.pos.sLoop()){
                v.x-=-a-b-d;
                v.y-=a-c-d;
                v.z-=b+c-d;
            }
        }
    } 
}

