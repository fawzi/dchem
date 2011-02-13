/// simplistic (convex) but very fast detection of fly away.
/// implement also neigh list/connectivity based analysis?
module dchem.pnet.DetectFlyAway;
import dchem.Common;
import dchem.pnet.PNetModels;
import dchem.input.RootInput;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.io.Console;
import blip.container.GrowableArray;
import dchem.sys.ParticleSys;
import dchem.input.WriteOut;
import dchem.pnet.EmptyObserver;
import dchem.sys.DynVars;
import dchem.sys.ParticleRange;
import blip.narray.NArray;
import blip.container.BulkArray;
import dchem.interpolate.Microsphere;
import dchem.util.Rotate;

class DetectFlyAwayGen:SilosWorkerGen{
    bool earlyDetect=true;
    Real threshold=3.0;
    InputField particles;
    this(){
    }
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.DetectFlyAway",`
    earlyDetect: if the configurations detected should be removed early (before trying to calculate them) (true)
    threshold: threshold that detects particles flying away (3.0)
    particles: the particles to monitor (if not given uses all the particles)`));
    mixin printOut!();
    bool verify(CharSink s){
        bool res=true;
        if (particles!is null && cast(ParticleRange)(particles.contentObj()) is null){
            dumper(s)("if particles is given it must be a ParticleRange in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
    
    SilosWorkerI!(Real) silosWorkerReal(){
        auto res=new DetectFlyAway!(Real)(this);
        return res;
    }
    SilosWorkerI!(LowP) silosWorkerLowP(){
        auto res=new DetectFlyAway!(LowP)(this);
        return res;
    }
}

/// an explorer that does nothing (useful as base class)
class DetectFlyAway(T):EmptyObserver!(T),SilosWorkerI!(T){
    DetectFlyAwayGen input;
    size_t lastNParticles=0;
    LocalSilosI!(T) silos;
    
    this(DetectFlyAwayGen input){
        this.input=input;
    }
    
    void workOn(LocalSilosI!(T) silos){
        silos=silos;
    }
    
    bool isFlyingAway(U=T)(ParticleSys!(T) pSys,DynPVector!(U,XType)xV){
        Vector!(U,3)[256] buf;
        auto xyz=lGrowableArray(buf,0,GASharing.Local);
        scope(exit) xyz.deallocData();
        if (input.particles!is null){
            foreach (pGroup;input.particles.contentT!(ParticleRange)().loopOn(pSys.sysStruct)){
                foreach(p;pGroup){
                    xyz.appendArr(xV.pos[p]);
                }
            }
        } else {
            foreach(k;xV.pos.kRange){
                xyz.appendArr(xV.pos[k].data);
            }
        }
        if (xyz.length==0) return false;
        auto xyzOrig=a2NA!(U,2)(bulkArray(xyz.data).basicData,false,true,[3,-1]);
        auto toSort=empty!(U)([3,xyz.length]);
        auto x=toSort[Range(0,-1),0].data;
        auto y=toSort[Range(0,-1),0].data;
        auto z=toSort[Range(0,-1),0].data;
        assert(x.length==xyz.length && y.length==xyz.length && z.length==xyz.length);
        auto dirs=default3SphereDirs.asType!(U)();
        auto startRot=default3SphereROrigin;
        foreach(idir,dir;dirs){
            toSort[]=xyzOrig;
            if (dir[0]!=1 || startRot[idir]!=0){
                rotateEiV(startRot[idir],dir,toSort);
            }
            sort(x);
            sort(y);
            sort(z);
            auto threshold=input.threshold;
            {
                auto oldV=x[0];
                for(size_t i=1;i<x.length;++i){
                    if (x[i]-oldV>threshold){
                        return true;
                    }
                }
            }
            {
                auto oldV=y[0];
                for(size_t i=1;i<y.length;++i){
                    if (y[i]-oldV>threshold){
                        return true;
                    }
                }
            }
            {
                auto oldV=z[0];
                for(size_t i=1;i<z.length;++i){
                    if (z[i]-oldV>threshold){
                        return true;
                    }
                }
            }
        }
        return false;
    }
    
    /// checks it local point is somehow invalid and should better be skipped
    bool shouldFilterLocalPoint(SKey s,Point p){
        if (input.earlyDetect){
            auto m=silos.mainPointL(p);
            alias isFlyingAway!(T) tt;
            if (isFlyingAway!(T)(m.pos,m.pos.dynVars.x)){
                return true;
            }
        }
        return false;
    }
    
    /// communicates that the given local point has been successfully published
    void publishedLocalPoint(SKey s,Point point){
        auto m=silos.mainPointL(point);
        if (isFlyingAway!(T)(m.pos,m.pos.dynVars.x)){
            m.gFlagsAtomicOp(delegate uint(uint a){ return (a|GFlags.DoNotExplore);});
        }
    }
    
}