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
import dchem.input.RootInput;
import dchem.sys.DynVars;

interface ConstraintGen:InputElement {
    ConstraintI!(Real) constraintReal(ParticleSys!(Real)templatePos);
    ConstraintI!(LowP) constraintLowP(ParticleSys!(LowP)templatePos);
}

/// represents one or more constraints
interface ConstraintI(T){
    /// the constraint generator that created this constraint
    ConstraintGen constraintGen();
    /// enforces the constraint on state, returns something 
    /// related to the change in the state (sup norm) + change in the contrained value
    /// if you have constraints on the velocities,... the update
    /// is not restricted to the x (position) values
    /// if partial is true only one iteration is performed, so the contraint might be unfulfilled
    /// (useful when called by a llop that tries to fulfull various constraints at once)
    real applyR(ParticleSys!(T) state,bool partial=false);
    /// removes the component of dr in state that is along a constraint, puts what is
    /// removed in diff if given. (projection state=(1-pi)state , diff=pi state, where pi
    /// is the projection in the space spanned by the gradient of the constraints).
    /// if you have non olonomic constraints or constraints on the velocities the update
    /// is not restricted to the dr (first derivatives) components
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null);
    /// value of a global constraint value of this form: sum((f_i(x)-f0_i)**2), where
    /// f_i are the various constraints, and f0_i the values they should be constrained to.
    /// will add the derivative of the constraint with respect to x to deriv.
    /// can be used to minimize the constraint error if iterating applyR (a la shake) has problems converging
    /// the value has to be positive, with 0 being the target
    real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null);
    /// iterates on the involvedparticles, there might be double counting, or extra particles
    /// on which the constraint is not really dependent
    FIteratorI!(PIndex)particlesInvolved();
    /// if the constraints are strictly on the positons of the particles listed
    bool strict();
}

/// helper to get a constraint of type T from a constraint generator
mixin(genTypeTMixin("Constraint","constraint","ParticleSys!(T) pSys","pSys"));

/// something that can loop on subsets of the particles
alias bool delegate(ref BulkArray!(PIndex)a) PSubsetLooper;

/// represents groups of particles.
/// The indexes of the particles (that are of level particleLevel) are referred with respect to a particle of kind
/// referencePKind in the level referenceLevel.
/// instead of the indexes it is possible to select all particles of a given kind (with particleKind).
/// The group repeats for all the particles of kind referencePKind. By default the reference is the whole system.
class ParticleRangeGen{
    PIndex[]particles;
    PIndex particleStart;
    PIndex particleEnd;
    long[] indexes;
    long indexStart=long.max;
    long indexEnd=long.min;
    char[] kindName;
    int kIdx=int.max;
    
    mixin(serializeSome("dchem.PRange",`
    particles: list of particles
    particleStart: first particle of a range that continues until particleEnd
    particleEnd: last particle of a range that starts at particleStart
    indexes: list of (external) indexes (starting at 0)
    indexStart: first (external) index of a range that ends at extIndexEnd
    indexEnd: last external (index) of a range that starts at extIndexStart
    kindName: name of the kind of particles that should be selected, if given alone selects all particles, otherwise sets the context for the other indexes (in that case it should be a molecule or chain name)
    kIdx: index of the kind of particles that should be selected, if given alone selects all particles, otherwise sets the context for the other indexes (in that case it should be a molecule or chain name)
    `));
    mixin printOut!();
    bool verify(CharSink sink,char[]fieldName){
        auto log=dumper(sink);
        bool res=true;
        bool hasParticles=false;
        bool hasIndex=false;
        if (particles.length!=0){
            hasParticles=true;
            if (particleStart.valid || particleEnd.valid) {
                log("either particles or particleStart and particleEnd have to be defined for ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
        }
        if (particleStart.valid || particleEnd.valid){
            hasParticles=true;
            if (!(particleStart.valid && particleEnd.valid)){
                log("Both particleStart and particleEnd or neither of them have to be defined for ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
        }
        if (indexes.length!=0){
            hasIndex=true;
            if (indexStart!=long.max || indexEnd!=long.min) {
                log("either indexes or indexStart and indexEnd have to be defined for ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
        }
        if (indexStart!=long.max || indexEnd!=long.min){
            hasIndex=true;
            if (!(indexStart!=long.max && indexEnd!=long.min)){
                log("Both indexStart and indexEnd or neither of them have to be defined for ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
        }
        if (hasIndex && hasParticles){
            log("indexes and particles cannot be defined together in ParticleRange in field ")
                (fieldName)("\n");
            res=false;
            
        }
        if (kindName.length!=0){
            if (kIdx!=int.max){
                log("kindName and kIdx cannot be defined together in ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
            if (hasParticles){
                log("if you define kindName only indexes can be used in ParticleRange in field ")
                    (fieldName)("\n");
                res=false;
            }
        }
        if (kIdx!=int.max && hasParticles){
            log("kIdx particles cannot be used (only indexes) in ParticleRange in field ")
                (fieldName)("\n");
            res=false;
        }
        return res;
    }
    
}

class ParticleRange(T){
    /+struct ElLooper(T){
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
    
    PSubsetLooper looperForPSys(ParticleSys!(T)pSys){
        auto l=loopOn(pSys);
        auto lNew=new typeof(l);
        *lNew=l;
        return &lNew.next;
    }+/
}

class MultiConstraintGen: ConstraintGen{
    InputField[] constraints;
    uint maxShakeIter;
    Real targetPrecision;

    mixin myFieldMixin!();
    mixin(serializeSome("dchem.MultiConstraint",`
    maxShakeIter: maximum iteration of "shake" (immediate partial inverse) before switching to direct optimization
    targetPrecision: the target precision that should be reached by the constraint loop
    constraints: the constraints to apply
    `));
    this(){}
    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        foreach (c;constraints){
            auto o=cast(Object)c.content;
            if (cast(ConstraintGen)o is null){
                log("The constraints have to be of type constraint, not ")((o !is null)?o.classinfo.toString():"*null*"[])(" in field ")(myField)("\n");
                res=false;
            }
        }
        return res;
    }
    ConstraintI!(Real) constraintReal(ParticleSys!(Real)templatePos){
        return constraintT!(Real)(templatePos);
    }
    ConstraintI!(LowP) constraintLowP(ParticleSys!(LowP)templatePos){
        return constraintT!(LowP)(templatePos);
    }
    
    MultiConstraint!(T) constraintT(T)(ParticleSys!(T)templatePos){
        return new MultiConstraint!(T)(this,templatePos);
    }
}

class MultiConstraint(T): ConstraintI!(T){
    MultiConstraintGen _constraintGen;
    ConstraintI!(T)[] subConstraints;
    
    this(MultiConstraintGen mc,ParticleSys!(T)templatePos){
        this._constraintGen=mc;
        subConstraints=new ConstraintI!(T)[](mc.constraints.length);
        foreach (i,c;mc.constraints){
            auto o=cast(ConstraintGen)cast(Object)c.content;
            assert(o!is null);
            subConstraints[i]=constraintT!(T)(o,templatePos);
        }
    }
    
    ConstraintGen constraintGen(){
        return _constraintGen;
    }
    /// builds overlap and checks if there are conflicting constraints
    bool hasConflicts(ParticleSys!(T) state,CharSink log){
        // to do build overlap matrix of the constraints, and check if there are conflicts
        return false;
    }
    
    real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null){
        real res=0;
        foreach(subC;subConstraints){
            res+=derivVal(state,deriv);
        }
        return res;
    }
    
    real applyR(ParticleSys!(T) state,bool partial=false){
        Real maxErr,oldMaxErr;
        bool checkConflicts=false;
        auto maxShakeIter=_constraintGen.maxShakeIter;
        auto targetPrecision=_constraintGen.targetPrecision;
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
    
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null){
        foreach(subC;subConstraints){
            subC.applyDR(state,diff);
        }
    }
    
    /// loops on all particles of all constraints
    static class ListParticles:FIteratorI!(PIndex){
        FIteratorI!(PIndex) iterAtt;
        ConstraintI!(T)[] c;
        size_t pos;
        this(ConstraintI!(T)[] c){
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

/// a constraint that applies no constraints
class NoConstraintGen:ConstraintGen {
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.NoConstraint",""));
    mixin printOut!();
    this(){}
    bool verify(CharSink s){
        return true;
    }
    ConstraintI!(Real) constraintReal(ParticleSys!(Real)templatePos){
        return new NoConstraint!(Real)(this);
    }
    ConstraintI!(LowP) constraintLowP(ParticleSys!(LowP)templatePos){
        return new NoConstraint!(LowP)(this);
    }
}
/// a constraint that applies no constraints
class NoConstraint(T):ConstraintI!(T){
    NoConstraintGen _constraintGen;
    this(){
        this._constraintGen=new NoConstraintGen();
    }
    this(NoConstraintGen c){
        this._constraintGen=c;
    }
    ConstraintGen constraintGen(){ return _constraintGen; }
    real applyR(ParticleSys!(T) state,bool partial=false){
        return 0;
    }
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null){ }
    real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null){ return 0; }
    FIteratorI!(PIndex)particlesInvolved(){ return EmptyFIterator!(PIndex).instance; }
    bool strict(){ return true; }
}

/// describes a function from some variables to another, the mapping should be derivable,
/// but not much more is required (i.e. it can be partial, surjective,...)
interface FunctionVarI(T){
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

interface FunctionVarGen{
    FunctionVarI!(Real) fVarReal(ParticleSys!(Real)templatePos);
    FunctionVarI!(LowP) fVarLowP(ParticleSys!(LowP)templatePos);
}

mixin(genTypeTMixin("FunctionVar","fVar","ParticleSys!(Real)templatePos","templatePos"));
/+
class DistanceConstraint{
}

class RigidBodyConstraint{
}

class FixedAtomConstraint{
}
+/
