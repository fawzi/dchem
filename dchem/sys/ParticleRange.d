/// describes ranges of particles 
module dchem.sys.ParticleRange;
import dchem.Common;
import blip.serialization.Serialization;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import blip.BasicModels;
import blip.container.BulkArray;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import dchem.input.RootInput;
import dchem.sys.DynVars;

/// an iterator that loops on subsets of the particles
alias bool delegate(ref BulkArray!(PIndex)a) PSubsetLooper;

/// represents groups of particles.
/// The indexes of the particles (that are of level particleLevel) are referred with respect to a particle of kind
/// referencePKind in the level referenceLevel.
/// instead of the indexes it is possible to select all particles of a given kind (with particleKind).
/// The group repeats for all the particles of kind referencePKind. By default the reference is the whole system.
class ParticleRange:InputField{
    PIndex[]particles;
    PIndex particleStart;
    PIndex particleEnd;
    long[] indexes;
    long indexStart=long.max;
    long indexEnd=long.min;
    char[] kindName;
    int kIdx=int.max;
    
    mixin(serializeSome("dchem.PRange",`
    particles: list of particles in internal notation (kind,)
    particleStart: first particle of a range that continues until particleEnd
    particleEnd: last particle of a range that starts at particleStart
    indexes: list of (external) indexes (starting at 0)
    indexStart: first (external) index of a range that ends at extIndexEnd
    indexEnd: last external (index) of a range that starts at extIndexStart
    kindName: name of the kind of particles that should be selected, if given alone selects all particles, otherwise sets the context for the other indexes (in that case it should be a molecule or chain name)
    kIdx: index of the kind of particles that should be selected, if given alone selects all particles, otherwise sets the context for the other indexes (in that case it should be a molecule or chain name)
    `));
    mixin printOut!();
    mixin myFieldMixin!();
    bool verify(CharSink sink){
        auto log=dumper(sink);
        bool res=true;
        bool hasParticles=false;
        bool hasIndex=false;
        if (particles.length!=0){
            hasParticles=true;
            if (particleStart.valid || particleEnd.valid) {
                log("either particles or particleStart and particleEnd but not both have to be defined for ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
        }
        if (particleStart.valid || particleEnd.valid){
            hasParticles=true;
            if (!(particleStart.valid && particleEnd.valid)){
                log("Both particleStart and particleEnd or neither of them have to be defined for ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
        }
        if (indexes.length!=0){
            hasIndex=true;
            if (indexStart!=long.max || indexEnd!=long.min) {
                log("either indexes or indexStart and indexEnd have to be defined for ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
        }
        if (indexStart!=long.max || indexEnd!=long.min){
            hasIndex=true;
            if (!(indexStart!=long.max && indexEnd!=long.min)){
                log("Both indexStart and indexEnd or neither of them have to be defined for ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
        }
        if (hasIndex && hasParticles){
            log("indexes and particles cannot be defined together in ParticleRange in field ")
                (myFieldName)("\n");
            res=false;
            
        }
        if (kindName.length!=0){
            if (kIdx!=int.max){
                log("kindName and kIdx cannot be defined together in ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
            if (hasParticles){
                log("if you define kindName only indexes can be used in ParticleRange in field ")
                    (myFieldName)("\n");
                res=false;
            }
        }
        if (kIdx!=int.max && hasParticles){
            log("kIdx particles cannot be used (only indexes) in ParticleRange in field ")
                (myFieldName)("\n");
            res=false;
        }
        return res;
    }
    
    static struct ElLooper(T){
        size_t posRefAtt;
        BulkArray!(PIndex) refParticles;
        ParticleSys!(T) pSys;
        long[] particleIndexes;
        BulkArray!(PIndex) realIndexes;
        /// visit the next reference particle (might invalidate the particleIndexes)
        bool next(ref BulkArray!(PIndex)a){
            if (posRefAtt<refParticles.length){
                if (particleIndexes){
                    auto subP=pSys.sysStruct.subParticles[refParticles[posRefAtt]];
                    foreach(i,v;particleIndexes){
                        realIndexes[i]=subP[cast(size_t)v];
                    }
                }
                a=realIndexes;
                ++posRefAtt;
                return true;
            }
            return false;
        }
        mixin opApplyFromNext!(BulkArray!(PIndex));
    }
    
    ElLooper!(T) loopOn(T)(ParticleSys!(T)pSys){
        ElLooper!(T) res;
        res.pSys=pSys;
        res.posRefAtt=0;
        if (referenceLevel==-1){
            res.refParticles=pSys.sysStruct.particles[pSys.sysStruct.levels[$-1].kStart];
            if (indexes.length!=0){
                auto extO=pSys.sysStruct.externalOrder;
                res.realIndexes=BulkArray!(PIndex)(indexes.length);
                foreach(i,p;indexes){
                    res.realIndexes[i]=PIndex(extO[LocalPIndex(cast(KindIdx)0,cast(ParticleIdx)p)]);
                }
            } else if (particleKind.length>0){
                bool found=false;
                foreach(pKind;pSys.sysStruct.particleKinds[pSys.sysStruct.levels[0]].sLoop){
                    if (pKind.name==particleKind){
                        found=true;
                        res.realIndexes=pSys.sysStruct.particles[pKind.pKind];
                        break;
                    }
                }
                if (!found){
                    throw new Exception(collectAppender(delegate void(CharSink sink){
                        auto s=dumper(sink);
                        s("could not find partikle kind named '")(particleKind)("', known particle kinds: ");
                        bool first=true;
                        foreach(pKind;pSys.sysStruct.particleKinds[pSys.sysStruct.levels[0]].sLoop){
                            if (first) s(", ");
                            first=false;
                            s("'")(pKind.name)("'");
                        }
                    }),__FILE__,__LINE__);
                }
            }
        } else {
            if (referenceLevel<1||referenceLevel>=pSys.sysStruct.levels.length){
                throw new Exception("referenceLevel out of bounds",__FILE__,__LINE__);
            }
            if (particlesLevel+1!=referenceLevel){
                throw new Exception("only system or references that are just one level above are implemented",__FILE__,__LINE__);
            }
            KindIdx kRef;
            foreach(pKind;pSys.sysStruct.particleKinds[pSys.sysStruct.levels[referenceLevel]].sLoop){
                if (pKind.name==referencePKind){
                    kRef=pKind.pKind;
                    break;
                }
            }
            if (kRef==KindIdx.init){
                throw new Exception(collectAppender(delegate void(CharSink sink){
                    auto s=dumper(sink);
                    s("could not find particle kind named '")(referencePKind)("', known particle kinds at level ")
                        (referenceLevel)(":");
                    bool nonFirst=false;
                    foreach(pKind;pSys.sysStruct.particleKinds[pSys.sysStruct.levels[referenceLevel]].sLoop){
                        if (nonFirst) s(", ");
                        nonFirst=true;
                        s("'")(pKind.name)("'");
                    }
                }),__FILE__,__LINE__);
            }
            res.refParticles=pSys.sysStruct.particles[kRef];
            
            if (particleIndexes.length!=0){
                res.particleIndexes=particleIndexes;
                res.realIndexes=BulkArray!(PIndex)(particleIndexes.length);
            } else if (particleKind.length>0){
                throw new Exception("particleKind supported only for global (system) references",__FILE__,__LINE__);
            }
        }

        return res;
    }
    
    PSubsetLooper looperForPSys(T)(ParticleSys!(T)pSys){
        auto lNew=new ElLooper!(T);
        *lNew=loopOn(pSys);
        return &lNew.next;
    }
}
