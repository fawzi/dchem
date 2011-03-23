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
import blip.io.FileStream;
import blip.core.Array;

/// detects fly away, at the moment it is convex, should I also implement the better connectivity based solution?
class DetectFlyAwayGen:ExplorationObserverGen{
    bool earlyDetect=false;
    char[] logBaseName;
    bool flushEachPoint=true;
    Real threshold=3.0;
    InputField particles;
    this(){
    }
    mixin myFieldMixin!();
    mixin(serializeSome("DetectFlyAway",`observer that detects when a system breaks apart and stop its exploration`,
    `earlyDetect: if the configurations detected should be removed early (before trying to calculate them) (false)
    logBaseName: if given logs the points that have "flown away" to a file that starts with logBaseName
    flushEachPoint: if each point should be immediately flushed (true)
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
    
    ExplorationObserverI!(Real) observerReal(LocalSilosI!(Real) silos){
        auto res=new DetectFlyAway!(Real)(this,silos);
        return res;
    }
    
    ExplorationObserverI!(LowP) observerLowP(LocalSilosI!(LowP) silos){
        auto res=new DetectFlyAway!(LowP)(this,silos);
        return res;
    }
}

/// an explorer that does nothing (useful as base class)
class DetectFlyAway(T):EmptyObserver!(T){
    DetectFlyAwayGen input;
    size_t lastNParticles=0;
    LocalSilosI!(T) silos;
    OutStreamI log;
    
    this(DetectFlyAwayGen input,LocalSilosI!(T) silos){
        this.input=input;
        this.silos=silos;
        if (input.logBaseName.length>0){
            log=outfileStrSync(input.logBaseName~"-"~silos.name~"-Filtred.Log",WriteMode.WriteAppend);
        }
    }
    
    bool isFlyingAway(U=T)(SKey sKey,Point point,ParticleSys!(T) pSys,DynPVector!(U,XType)xV){
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
        U[][3] xyzSplit;
        for (int idim=0;idim<3;++idim){
            xyzSplit[idim]=toSort[idim].data;
            assert(xyzSplit[idim].length==xyz.length);
        }
        auto dirs=default3SphereDirs.asType!(U)();
        auto startRot=default3SphereROrigin;
        foreach(idir,dir;dirs){
            toSort[]=xyzOrig;
            if (dir[0]!=1 || startRot[idir]!=0){
                rotateEiV(startRot[idir],dir,toSort);
            }
            auto threshold=input.threshold;
            for (int idim=0;idim<3;++idim){
                sort(xyzSplit[idim]);
                auto xx=xyzSplit[idim];
                auto oldV=xx[0];
                for(size_t i=1;i<xx.length;++i){
                    if (xx[i]-oldV>threshold){
                        if (log!is null){
                            auto mDir=zeros!(U)(3);
                            mDir[idim]=cast(U)1;
                            rotateEiV(startRot[idir],dir,mDir);
                            auto mid=(xx[i]+oldV)/2;
                            sinkTogether(log.charSink(),delegate void(CharSink s){
                                dumper(s)(sKey)(" \t")(point)(" \t")(mDir[0])(" \t")(mDir[1])(" \t")(mDir[2])(" \t")
                                    (mid*mDir[0])(" \t")(mid*mDir[1])(" \t")(mid*mDir[2])("\n");
                            });
                        }
                        return true;
                    }
                    oldV=xx[i];
                }
            }
        }
        return false;
    }
    
    /// checks it local point is somehow invalid and should better be skipped
    override bool shouldFilterLocalPoint(SKey s,Point p){
        if (input.earlyDetect){
            auto m=silos.mainPointL(p);
            alias isFlyingAway!(T) tt;
            if (isFlyingAway!(T)(s,p,m.pos,m.pos.dynVars.x)){
                return true;
            }
        }
        return false;
    }
    
    /// communicates that the given local point has been successfully published
    override void publishedLocalPoint(SKey s,Point point){
        auto m=silos.mainPointL(point);
        if (isFlyingAway!(T)(s,point,m.pos,m.pos.dynVars.x)){
            m.gFlagsAtomicOp(delegate uint(uint a){ return (a|GFlags.DoNotExplore);});
        }
    }
    
    override void increaseRunLevel(SKey s,RunLevel level){
        if (log!is null){
            log.flush();
            if (level>=RunLevel.Stopped){
                log.close();
            }
        }
    }
}
