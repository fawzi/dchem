module dchem.sys_struct.AtomParticle;

/// represent atom like particleKindKind (has no direct pointer to actual particles)
class AtomicKindKind: ParticleKindKind,CopiableOjectI,Serializable {
    char[] _particleKindName;
    char[] atomSymbol;
    PropertyKind mass,pos,velocity,force,charge;
    PropertyKind[] _propertyKinds;
    char[] particleKindName() {
        return _particleKindName;
    }
    int level(){
        return 0;
    }
    bool fixedLayout(){
        true;
    }
    PropertyKind[] propertyKinds(){
        _propertyKinds;
    }
    PropertyKind getPropertyKind(char[] propertyName){
        switch(propertyName){
            case "mass": return m;
            case "pos": return pos;
            case "velocity": return velocity;
            case "force": return force;
            case "charge": return charge;
            default:
                foreach(p;_propertyKinds){
                    if (p.name==propertyName) return p;
                }
                return null;
            
        }
    }
    void serialize(Serializer s);
    ParticleKindKind dup(PSCopyDepthLevel level);
}


/// represent a group of particles with the same kind
interface ParticleKind: CopiableOjectI,Serializable{
    char[] particleKindName();
    int level();
    int pKindIdx();
    GenProperty[] properties();
    GenProperty propertyNamed(char[] propertyName);

    BulkArray!(idxType) particlesArray();
    void particlesArray(BulkArray!(idxType) pArray);
    FIteratorI!(idxType) particlesIdxs();
    FIteratorI!(Particle *) particles();
    
    // sub particles
    FIteratorI!(idxType) subParticleIdxs(idxType localParticleIdx);
    idxType subParticleIdx(idxType localParticleIdx,idxType subParticleIdx);
    FIteratorI!(Particle*) subParticles(idxType localParticleIdx);
    Particle* subParticle(idxType localParticleIdx,idxType subParticleIdx);
    
    void serialize(Serializer s);
    ParticleKindKind kind();
    bool generateSubparticleIndexes(); /// if the level info should keep indexes about subParticles
    ParticleKindKind dup(PSCopyDepthLevel level);
}
