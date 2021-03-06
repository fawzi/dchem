/// data stored for each local point, storage efficient for data that needs to be stored for each point (dense)
/// useful to unbundle extra method data from the points
module dchem.pnet.DenseLocalPointArray;
import blip.container.BatchedGrowableArray;
import blip.container.UnserializableBatchedGrowableArray;
import blip.serialization.Serialization;
import dchem.pnet.PNetModels;
import blip.util.BinSearch;
import blip.sync.Atomic;

/// dense array of elements indexed by local point keys
class DenseLocalPointArray(T){
    BatchedGrowableArray!(Point,batchSize) keys;
    UnserializableBatchedGrowableArray!(T,batchSize) values;
    
    T opIndex(Point k){
        synchronized(this){
            auto idx=lBound(keys,k);
            assert(idx<keys.length && keys[idx]==k,"could not find local key "~k.toString());
            if (values.length<idx){
                return T.init;
            }
            return values[idx];
        }
    }
    T* ptrI(Point k){
        size_t nKeys=keys.length;
        readBarrier();
        auto idx=lBound(keys,k,0,nKeys);
        assert(idx<nKeys && keys[idx]==k,"could not find local key "~k.toString());
        if (values.length<=idx){
            values.growTo(keys.length);
        }
        auto res=values.ptrI(idx);
        return res;
    }
    void opIndexAssign(T val,Point k){
        synchronized(this){
            auto vPos=ptrI(k);
            *vPos=val;
        }
    }
    T* opIn_r(Point k){
        size_t nKeys=keys.length;
        readBarrier();
        auto idx=lBound(keys,k,0,nKeys);
        if (idx<nKeys && keys[idx]==k){
            if (values.length<idx){
                values.growTo(keys.length);
            }
            return values.ptrI(idx);
        }
        return null;
    }
    
    int opApply(int delegate(ref Point,ref T el)loopBody){
        auto ks=keys.view();
        synchronized(this){
            values.growTo(ks.length);
        }
        foreach(i,k;ks){
            auto res=loopBody(k,*values.ptrI(i));
            if (res!=0) return res;
        }
        return 0;
    }
    
    size_t length(){
        return keys.length;
    }
    
    this(BatchedGrowableArray!(Point,batchSize) pointsKey){
        keys=pointsKey;
        values=new UnserializableBatchedGrowableArray!(T,batchSize)();
    }
    this(){
        this(new BatchedGrowableArray!(Point,batchSize)());
    }
}
