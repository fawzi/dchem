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
import dchem.sys.ParticleRange;

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
    Real applyR(ParticleSys!(T) state,bool partial=false);
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
    Real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null);
/+    /// iterates on the involvedparticles, there might be double counting, or extra particles
    /// on which the constraint is not really dependent
    FIteratorI!(PIndex)particlesInvolved();
    /// if the constraints are strictly on the positons of the particles listed
    bool strict();+/
}

/// helper to get a constraint of type T from a constraint generator
mixin(genTypeTMixin("Constraint","constraint","ParticleSys!(T) pSys","pSys"));

class MultiConstraintGen: ConstraintGen{
    InputField[] constraints;
    uint maxShakeIter;
    Real targetPrecision=1.0e-6;

    mixin myFieldMixin!();
    mixin(serializeSome("dchem.MultiConstraint",`Constraint that Handles a group of constraints.`,
    `maxShakeIter: maximum iteration of "shake" (immediate partial inverse) before switching to direct optimization
    targetPrecision: the target precision that should be reached by the constraint loop
    constraints: the constraints to apply
    `));
    mixin printOut!();
    
    this(){}
    bool verify(CharSink s){
        auto log=dumper(s);
        bool res=true;
        foreach (c;constraints){
            auto o=cast(Object)c.content;
            if (cast(ConstraintGen)o is null){
                log("The constraints have to be of type constraint, not ")((o !is null)?o.classinfo.name:"*null*"[])(" in field ")(myFieldName)("\n");
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
    
    Real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null){
        Real res=0;
        foreach(subC;subConstraints){
            res+=derivVal(state,deriv);
        }
        return res;
    }
    
    Real applyR(ParticleSys!(T) state,bool partial=false){
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
        return maxErr;
    }
    
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null){
        foreach(subC;subConstraints){
            subC.applyDR(state,diff);
        }
    }
    
/+    /// loops on all particles of all constraints
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
    }+/
}

/// a constraint that applies no constraints
class NoConstraintGen:ConstraintGen {
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.NoConstraint","Empty constraint.",""));
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
    Real applyR(ParticleSys!(T) state,bool partial=false){
        return 0;
    }
    void applyDR(ParticleSys!(T) state,ParticleSys!(T) diff=null){ }
    Real derivVal(ParticleSys!(T) state,ParticleSys!(T) deriv=null){ return 0; }
/+    FIteratorI!(PIndex)particlesInvolved(){ return EmptyFIterator!(PIndex).instance; }
    bool strict(){ return true; }+/
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

/+ // to do
class DistanceConstraint{
}

class RigidBodyConstraint{
}

class FixedAtomConstraint{
}
+/
