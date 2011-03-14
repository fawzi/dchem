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
class ParticleRange:InputElement{
    PIndex[]particles;
    PIndex particleStart;
    PIndex particleEnd;
    long[] indexes;
    long indexStart=long.max;
    long indexEnd=long.min;
    char[] kindName;
    int kIdx=int.max;
    long refKind=-1;
    char[] refKindName;
    
    mixin(serializeSome("dchem.PRange",`Structure that represent some particles (subsets,...)`,
    `particles: list of particles in internal notation (kind,idx)
    particleStart: first particle of a range that continues until particleEnd
    particleEnd: last particle of a range that starts at particleStart
    indexes: list of (external) indexes (starting at 0)
    indexStart: first (external) index of a range that ends at extIndexEnd
    indexEnd: last (external) index of a range that starts at extIndexStart
    kIdx: index of the kind of particles that should be selected, if given alone selects all particles, otherwise sets the context for the other indexes
    kindName: name of the kind of particles that should be selected, if given alone selects all particles
    refKind: the level of the reference (by default -1, i.e. system level)
    refKindName: the name of the level of the reference (overrides refKind)
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
    
    static struct ElLooper{
        size_t posRefAtt;
        BulkArray!(PIndex) refParticles;
        SysStruct sysStruct;
        long[] particleIndexes;
        BulkArray!(PIndex) realIndexes;
        /// visit the next reference particle (might invalidate the particleIndexes)
        bool next(ref BulkArray!(PIndex)a){
            if (posRefAtt<refParticles.length){
                if (particleIndexes){
                    auto subP=sysStruct.subParticles[refParticles[posRefAtt]];
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
    
    ElLooper loopOn(SysStruct sysStruct){
        ElLooper res;
        res.sysStruct=sysStruct;
        res.posRefAtt=0;
        auto referenceK=refKind;
        if (refKindName.length>0){
            bool found=false;
            foreach(pKind;sysStruct.particleKinds[sysStruct.levels[0]].sLoop){
                if (pKind.name==refKindName){
                    found=true;
                    referenceK=pKind.pKind;
                    break;
                }
            }
            if (!found) throw new Exception("could not find refKindName '"~refKindName~"' in field "~myFieldName,
                __FILE__,__LINE__);
        }
        if (referenceK==-1 || referenceK==sysStruct.levels.length-1){
            res.refParticles=sysStruct.particles[sysStruct.levels[$-1].kStart];
            if (indexes.length!=0|| indexStart<indexEnd){
                auto extO=sysStruct.externalOrder;
                auto len=indexes.length;
                if (indexStart<indexEnd){
                    len+=cast(size_t)indexEnd-indexStart;
                }
                res.realIndexes=BulkArray!(PIndex)(len);
                foreach(i,p;indexes){
                    res.realIndexes[i]=PIndex(extO[PIndex(cast(KindIdx)0,cast(ParticleIdx)p)]);
                }
                size_t ii=indexes.length;
                for (long i=indexStart;i<indexEnd;++i){
                    res.realIndexes[ii++]=PIndex(extO[PIndex(cast(KindIdx)0,cast(ParticleIdx)i)]);
                }
                assert(kindName.length==0);
                if (kindName.length>0){
                    throw new Exception("kindName supported only without indexes",__FILE__,__LINE__);
                }
            } else if (kindName.length>0){
                bool found=false;
                foreach(pKind;sysStruct.particleKinds[sysStruct.levels[0]].sLoop){
                    if (pKind.name==kindName){
                        found=true;
                        res.realIndexes=sysStruct.particles[pKind.pKind];
                        break;
                    }
                }
                if (!found){
                    throw new Exception(collectAppender(delegate void(CharSink sink){
                        auto s=dumper(sink);
                        s("could not find partikle kind named '")(kindName)("', known particle kinds: ");
                        bool first=true;
                        foreach(pKind;sysStruct.particleKinds[sysStruct.levels[0]].sLoop){
                            if (first) s(", ");
                            first=false;
                            s("'")(pKind.name)("'");
                        }
                    }),__FILE__,__LINE__);
                }
            }
        } else {
            if (referenceK<1||referenceK>=sysStruct.levels.length){
                throw new Exception("referenceK out of bounds",__FILE__,__LINE__);
            }
            res.refParticles=sysStruct.particles[cast(KindIdx)referenceK];
            if (indexes.length>0){
                res.particleIndexes=indexes;
            }
            if (indexStart<indexEnd){
                auto nRange=cast(size_t)(indexEnd-indexStart);
                res.particleIndexes.length=res.particleIndexes.length+nRange;
                foreach (i,ref v;res.particleIndexes[res.particleIndexes.length-nRange..$]){
                    v=indexStart+i;
                }
            }
            res.realIndexes=BulkArray!(PIndex)(res.particleIndexes.length);
            if (kindName.length>0){
                throw new Exception("kindName supported only for global (system) references",__FILE__,__LINE__);
            }
        }
        return res;
    }
    
    PSubsetLooper looperForSysStruct(SysStruct sysStruct){
        auto lNew=new ElLooper;
        *lNew=loopOn(sysStruct);
        return &lNew.next;
    }
}
