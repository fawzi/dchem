/// data stored for each local point, storage efficient for data that needs to be stored for each point (dense)
/// useful to unbundle extra method data from the points
module dchem.pnet.DenseLocalPointArray;

enum { batchSize=512 }

/// dense array of elements indexed by local point keys
class DenseLocalPointArray(T){
    BatchedGrowableArray!(Point,batchSize) keys;
    BatchedGrowableArray!(T,batchSize) values;
    
    T opIndex(Point k){
        synchronized(this){
            auto idx=lbound(keys,k);
            assert(idx<keys.length && keys[idx]==k,"could not find local key "~k.toString());
            if (values.length<idx){
                return T.init;
            }
            return values[idx];
        }
    }
    T* ptrI(Point k){
        auto idx=lbound(keys,k);
        assert(idx<keys.length && keys[idx]==k,"could not find local key "~k.toString());
        if (values.length<idx){
            values.growTo(keys.length);
        }
        return &(values[idx]);
    }
    void opIndexAssign(T val,Point k){
        synchronized(this){
            *ptrI()=val;
        }
    }
    /+static if (T U:U*){
        int opApply(int delegate(ref U el)){
            synchronized(this){
                auto k=keys.view();
                values.growTo(k.length);
                auto v=values.view();
            }
            //v.sLoop()
        }
    }+/
}
