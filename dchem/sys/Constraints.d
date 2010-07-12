module dchem.sys.Constraints;
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

/// represents one or more constraints
interface Constraint{
    /// enforces the constraint on state, returns something 
    /// related to the change in the state (sup norm) + change in the contrained value
    /// if you have constraints on the velocities,... the update
    /// is not restricted to the x (position) values
    /// if partial is true only one iteration is performed, so the contraint might be unfulfilled
    /// (useful when called by a llop that tries to fulfull various constraints at once)
    real applyR(ParticleSys!(Real) state,bool partial=false);
    /// ditto
    real applyR(ParticleSys!(LowP) state,bool partial=false);
    /// removes the component of dr in state that is along a constraint, puts what is
    /// removed in diff if given. (projection state=(1-pi)state , diff=pi state, where pi
    /// is the projection in the space spanned by the gradient of the constraints).
    /// if you have non olonomic constraints or constraints on the velocities the update
    /// is not restricted to the dr (first derivatives) components
    void applyDR(ParticleSys!(Real) state,ParticleSys!(Real) diff=null);
    /// ditto
    void applyDR(ParticleSys!(LowP) state,ParticleSys!(LowP) diff=null);
    /// value of a global constraint value of this form: sum((f_i(x)-f0_i)**2), where
    /// f_i are the various constraints, and f0_i the values they should be constrained to.
    /// will add the derivative of the constraint with respect to x to deriv.
    /// can be used to minimize the constraint error if iterating applyR (a la shake) has problems converging
    /// the value has to be positive, with 0 being the target
    real derivVal(ParticleSys!(Real) state,ParticleSys!(Real) deriv=null);
    /// ditto
    real derivVal(ParticleSys!(LowP) state,ParticleSys!(LowP) deriv=null);
    /// iterates on the involvedparticles, there might be double counting, or extra particles
    /// on which the constraint is not really dependent
    FIteratorI!(PIndex)particlesInvolved();
    /// if the constraints are strictly on the positons of the particles listed
    bool strict();
}

/// something that can loop on subsets of the particles
alias bool delegate(ref BulkArray!(PIndex)a) PSubsetLooper;

interface ParticleSelector:Serializable{
    PSubsetLooper looperForPSys(ParticleSys!(Real));
    PSubsetLooper looperForPSys(ParticleSys!(LowP));
}

/// represents groups of particles.
/// The indexes of the particles (that are of level particleLevel) are referred with respect to a particle of kind
/// referencePKind in the level referenceLevel.
/// instead of the indexes it is possible to select all particles of a given kind (with particleKind).
/// The group repeats for all the particles of kind referencePKind. By default the reference is the whole system.
struct Particles{
    int referenceLevel=-1;
    char[] referencePKind="System";
    int particlesLevel=0;
    ulong[] particleIndexes;
    char[] particleKind;
    
    mixin(serializeSome("dchem.Particles",
    `referenceLevel: the level of the reference particle (defaults to -1, the whole system)
    referencePKind: the name of the reference particleKind
    particleIndexes: indexes of the particles
    particleKind: instead of the indexes the kind of the particles can be specified`));
    mixin printOut!();
    bool verify(CharSink sink,char[]fieldName){
        if (particleKind.length!=0 && particleIndexes.length!=0){
            dumper(sink)("Error: only one of particleIndexes and particleKind should be set in field "~fieldName);
            return false;
        }
    }
    
    struct ElLooper(T){
        size_t posRefAtt;
        BulkArray!(PIndex) refParticles;
        ParticleSys!(T) pSys;
        ulong[] particleIndexes;
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
    }
    
    ElLooper!(T) loopOn(T)(ParticleSys!(T)pSys){
        ElLooper!(T) res;
        res.pSys=pSys;
        res.posRefAtt=0;
        if (referenceLevel==-1){
            res.refParticles=pSys.sysStruct.particles[pSys.sysStruct.levels[$-1].kStart];
            
            if (particleIndexes.length!=0){
                auto extO=pSys.sysStruct.externalOrder;
                res.realIndexes=BulkArray!(PIndex)(particleIndexes.length);
                foreach(i,p;particleIndexes){
                    res.realIndexes[i]=extO[LocalPIndex(p)];
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
    
    PSubsetLooper looperForPSys(ParticleSys!(Real)pSys){
        auto l=loopOn(pSys);
        auto lNew=new typeof(l);
        *lNew=l;
        return &lNew.next;
    }
    PSubsetLooper looperForPSys(ParticleSys!(LowP)pSys){
        auto l=loopOn(pSys);
        auto lNew=new typeof(l);
        *lNew=l;
        return &lNew.next;
    }
}

class MultiConstraint: Constraint{
    Constraint[] subConstraints;
    uint maxShakeIter;
    Real targetPrecision;
    
    /// builds overlap and checks if there are conflicting constraints
    bool hasConflicts(T)(ParticleSys!(T) state,CharSink log){
        // to do build overlap matrix of the constraints, and check if there are conflicts
        return false;
    }
    
    real derivValT(T)(ParticleSys!(T) state,ParticleSys!(T) deriv=null){
        real res=0;
        foreach(subC;subConstraints){
            res+=derivVal(state,deriv);
        }
        return res;
    }
    
    real applyRT(T)(ParticleSys!(T) state,bool partial=false){
        Real maxErr,oldMaxErr;
        bool checkConflicts=false;
        for (uint i=0;i!=maxShakeIter;++i){
            maxErr=0;
            foreach (subC;subConstraints){
                auto err=subC.applyR(state,true);
                maxErr=max(maxErr,err);
            }
            if (maxErr<targetPrecision || partial) break;
            if (i!=0) {
                if (maxErr>=oldMaxErr){
                    checkConflicts=true;
                    break;
                }
            }
            oldMaxErr=maxErr;
        }
        if (checkConflicts){
            hasConflicts(state,serr.call); // should probably use another log...
        }
        if (!partial && maxErr>targetPrecision){
            // to do cg optimization
            assert(0,"cg optimization not yet implemented");
        }
    }
    real applyR(ParticleSys!(LowP) state,bool partial=false){
        return applyRT!(LowP)(state,partial);
    }
    real applyR(ParticleSys!(Real) state,bool partial=false){
        return applyRT!(Real)(state,partial);
    }
    
    void applyDRT(T)(ParticleSys!(T) state,ParticleSys!(T) diff=null){
        foreach(subC;subConstraints){
            subC.applyDR(state,diff);
        }
    }
    void applyDR(ParticleSys!(Real) state,ParticleSys!(Real) diff=null){
        applyDRT!(Real)(state,diff);
    }
    void applyDR(ParticleSys!(LowP) state,ParticleSys!(LowP) diff=null){
        applyDR!(LowP)(state,diff);
    }
    
    real derivVal(ParticleSys!(Real) state,ParticleSys!(Real) deriv=null){
        return derivValT(state,deriv);
    }

    real derivVal(ParticleSys!(LowP) state,ParticleSys!(LowP) deriv=null){
        return derivValT(state,deriv);
    }
    
    /// loops on all particles of all constraints
    static class ListParticles:FIteratorI!(PIndex){
        FIteratorI!(PIndex) iterAtt;
        Constraint[] c;
        size_t pos;
        this(Constraint[] c){
            this.c=c;
            pos=0;
            iterAtt=null;
        }
        bool next(ref PIndex t){
            while(true){
                if(iterAtt !is null){
                    if (iterAtt.next(t)){
                        return true;
                    }
                }
                if (pos<c.length){
                    iterAtt=c[pos].particlesInvolved();
                    ++pos;
                } else {
                    return false;
                }
            }
        }
        mixin opApplyFromNext!(PIndex);
        ForeachableI!(PIndex) parallelLoop(size_t optimalChunkSize){
            return this;
        }
        ForeachableI!(PIndex) parallelLoop(){
            return this;
        }
    }
    /// iterates on the involvedparticles, there might be double counting
    FIteratorI!(PIndex)particlesInvolved(){
        auto res=new ListParticles(subConstraints);
        return res;
    }

    /// if the constraints are strictly on the positons of the particles listed
    bool strict(){
        bool res=true;
        foreach (subC;subConstraints){
            res=res&&subC.strict();
        }
        return res;
    }
}

/// describes a function from some variables to another, the mapping should be derivable,
/// but not much more is required (i.e. it can be partial, surjective,...)
/// (this should be a templated interface that thing would be less buggy)

// f()PSys->DynPVect2 npos,norient,ndof
// particles
interface FunctionVar(T){
    void addValueTo(ParticleSys!(T) pos,DynPVector!(T,XType)target);
    /// transfers the derivative functionDeriv with respect to this variables to derivatives with respect
    /// to the pos variable, and adds them to targetTs:
    /// targetTs+=functionDeriv*df/dpos|_{pos}
    void addDerivToTs(T)(ParticleSys!(T) pos,DynPVector!(T,DxType)functionDeriv,DynPVector!(T,DxType)targetTs);
    /// performs the scalar product of the vector p with the derivative of this function
    /// i.e. it returns dot(functionDeriv*df/dpos|_{pos},vect)
    T collectDerivInTs(T)(ParticleSys!(T) pos,DynPVector!(T,DxType)functionDeriv,DynPVector!(T,DualDxType) vect);
    /// tries to set pos, so that evaluating it one obtains fromT
    /// returns true if it was successful, false if it was only partially done, or not at all
    bool quickInverse(T)(ParticleSys!(T) pos,DynPVector!(T,XType) fromT);
}
/+
class DistanceConstraint{
}

class RigidBodyConstraint{
}

class FixedAtomConstraint{
}
+/
