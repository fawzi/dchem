module blip.chem.sysStruct.PProperties;
import blip.narray.NArray;
import blip.Serializer

/// a particle property using an NArray
class NArrayPProperty!(T,int rank,bool readOnly=false): ParticleProperty!(NArray!(T,rank)){
    /// creates a property storage with the given shape and for the given number of particles
    this(PropertyKind pKind,index_type[rank] shape,idxType nParticles){
        index_type[rank+1] shp;
        shp[0]=nParticles;
        shp[1..$]=shape;
        this.arr=zeros!(T)(shp);
        this._pKind=pKind;
    }
    PropertyKind _pKind;
    PropertyKind pKind(){ return _pKind; }
    /// the type of the current property
    alias T PType;
    /// the array storing the property
    NArray!(T,rank+1) arr;
    /// get the value for the particle at index i
    NArray!(T,rank) opIndex(idxType i){
        return arr.opIndex(i);
    }
    /// sets the value for the particle at index i
    void opIndexAssign(idxType i,NArray!(T,rank) val){
        arr.opIndexAssign(i,val);
    }
    /// sets the value for all the particles
    void value(NArray!(T,rank) val){
        arr[]=repeat!(T,rank)(val,arr.shape[0]);
    }
}

/// a particleKind property using an NArray
class NArrayPKindProperty!(T,int rank,bool readOnly=false): PKindProperty!(NArray!(T,rank)){
    /// creates a property storage with the given shape and for the given number of particles
    this(PropertyKind pKind,index_type[rank] shape,idxType nParticles){
        this.arr=zeros!(T)(shape);
        this._pKind=pKind;
    }
    PropertyKind _pKind;
    PropertyKind pKind(){ return _pKind; }
    /// the type of the current property
    alias T PType;
    /// the array storing the property
    NArray!(T,rank) arr;
    /// get the value for the particle at index i
    NArray!(T,rank) opIndex(idxType i){
        return arr
    }
    /// sets the value for all the particles
    void value(NArray!(T,rank) val){
        arr[]=val;
    }
    NArray!(T,rank) value(){
        return arr;
    }
}
