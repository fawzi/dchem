module dchem.parallel.Distribution1d;

class Distribution1d: Serializable{
    /// communicator
    LinearComm comm;
    /// owner of particle i
    SegmentedArray!(int) owner;
    /// submap to local particles
    SubMapping localParticles;
    
    this(){ assert(false,"creation should normally at least assign a linear communicator"); }
    this(LinearComm comm){
        
    }
    this(LinearComm comm, SegmentedArray!(int)owner,SubMapping localParticles){
        
    }
    
}

