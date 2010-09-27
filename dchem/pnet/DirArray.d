/// storage of flags to keep track of explored directions.
/// this did change during development, as it went to single bit. It should be rewritten using BitArray
/// (that I already began adapting, but for now it still uses the tested code based on 2 bit per direction, 
/// as it was well tested.
module dchem.pnet.DirArray;
import blip.rtest.BasicGenerators;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.container.GrowableArray;
import Atomic=blip.sync.Atomic;

/// flags about a direction
/// after the rewrite there is one bit free here...
enum DirFlags{
    Free=0, /// the direction has not been explored
    Explored=1, /// the direction is being explored
    BlockedByMethod=2, /// declared as explored because the current method does not want to explore this direction
    BlockedStruct=3 /// declared as explored because constraints, coordinates choice or potential energy surface makes us skip it
}

/// prints the bit pattern of a, at least minBits are printed (minBits can be at most 64)
char[] printBits(size_t a,int minBits=1){
    char[64] res;
    auto i=64;
    while(a!=0){
        --i;
        if ((a&1)==1){
            res[i]='1';
        } else {
            res[i]='0';
        }
        a>>=1;
    }
    auto mBit=((minBits<=64)?minBits:64);
    while (64-i<mBit){
        res[--i]='0';
    }
    return res[i..$].dup;
}

/// an array with flags for each direction, does not use BitArray to have a fast earch of free directions
/// in large almost full arrays
///
/// At the moment it seems that only one bit per direction would be enough with the current setup, but
/// I am waiting for removing the support for two bits as it is tested and not fully trivial...
class FlagsArray{
    enum{ 
        usedBits=2,
        bitsPerEl=size_t.sizeof*8/usedBits,
        mask=((~(cast(size_t)0))>>(size_t.sizeof*8-usedBits))
    }
    static if (size_t.sizeof==4){
        enum :size_t{
            add0        =0b10101010_10101010_10101010_10101011,
            bitToCheck  =0b10101010_10101010_10101010_10101010,
        }
    } else {
        enum :size_t{
            add0        =0b10101010_10101010_10101010_10101010_10101010_10101010_10101010_10101011,
            bitToCheck  =0b10101010_10101010_10101010_10101010_10101010_10101010_10101010_10101010,
        }
    }
    size_t nFlags;
    size_t nCheckLast;
    size_t[] data;
    this(){}
    /// construct an array with the given size
    this(size_t dim){
        nFlags=dim;
        data=new size_t[]((dim+bitsPerEl-1)/bitsPerEl);
    }
    /// generates a random FlagsArray array
    static FlagsArray randomGenerate(Rand r,ref bool acceptable){
        auto dim=generateSize(r,size_t.sizeof);
        auto res=new FlagsArray(dim);
        mkRandomArray(r,res.data,acceptable);
        return res;
    }
    // returns the flags at the given index
    uint opIndex(size_t idx)
    in{
        if (idx>=nFlags){
            throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)("flags index ")(idx)(" out of bounds (0,")(nFlags)(")");
            }),__FILE__,__LINE__);
        }
    } body {
        auto block=data[idx/bitsPerEl];
        block>>=usedBits*(idx%bitsPerEl);
        return (mask&block);
    }
    /// atomic cas on a flag 
    uint atomicCAS(size_t idx,uint newF,uint oldF){
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        auto oldVal=data[blockIdx];
        assert((newF&mask)==newF,"new value is too large");
        auto localVal=(cast(size_t)newF)<<inBlockIdx;
        while(1){ // only for collisions with other updates
            auto oldFlags=((oldVal>>inBlockIdx)&0b11);
            if (oldFlags!=oldF) return oldFlags;
            auto actualVal=Atomic.atomicCAS(data[blockIdx],((oldVal&inBlockMask)|localVal),oldVal);
            if ((actualVal&inBlockMask)==(oldVal&inBlockMask)) return ((actualVal>>inBlockIdx)&0b11);
            oldVal=actualVal;
        }
    }
    /// atomic operation on a flag 
    void atomicOp(size_t idx,uint delegate(uint) op){
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        Atomic.atomicOp(data[blockIdx],delegate size_t(size_t oldVal){
            auto oldFlags=(oldVal>>inBlockIdx)&0b11;
            auto val=op(oldVal);
            assert((val&mask)==val,"val is too large");
            auto localVal=(cast(size_t)val)<<inBlockIdx;
            return ((oldVal&inBlockMask)|localVal);
        });
    }
    /// sets the flags at the given index (atomic, updates to *other* indexes cannot clash)
    void opIndexAssign(uint val,size_t idx){
        assert((val&mask)==val,"val is too large");
        assert(idx<nFlags,"flags index out of bounds");
        auto blockIdx=idx/bitsPerEl;
        auto inBlockIdx=(idx%bitsPerEl)*usedBits;
        auto inBlockMask=~(mask<<inBlockIdx);
        auto localVal=((cast(size_t)val)<<inBlockIdx);
        atomicOp(data[blockIdx],delegate size_t(size_t oldVal){
            return ((oldVal&inBlockMask)|localVal);
        });
    }
    /// returns a value different from 0 if the value val contains 00 (aligned)
    /// 4 ops + compare to 0 -> checks that 16 or 32 positions are not 00
    /// some nice bit trickery by me (fawzi :), actually the positions of the 1 bit corresponds to
    /// the start of the 00 in all cases but xxx0100 that is mapped to yyy1010 
    /// (but I did not try to take advantage of that)
    static size_t contains00(size_t val){
        auto r=((((val+add0)|val)^val)&bitToCheck);
        return r;
    }
    /// looks for a free (00) index starting at startPos.
    /// wraps around, looks at nFlags-1 as last direction, returns nFlags if no direction was found
    /// if last last is true then the last direction is checked in any case as last
    size_t findFreeAndSet(size_t startPos,uint newVal=0,bool lastLast=true,FlagsArray extraCheck=null){
        if (nFlags==0) return 0;
        if (nFlags==1) {
            if (atomicCAS(0,newVal,0)==0){
                return 0;
            }
            return 1;
        }
        assert(startPos<nFlags,"startPos out of bounds");
        assert(extraCheck is null || extraCheck.length>=length,"extraCheck is too short");
        auto idx=startPos;
        auto blockIdx=idx/bitsPerEl;
        auto rest=idx%bitsPerEl;
        auto val=(data[blockIdx])>>(rest*usedBits);
        rest=bitsPerEl-rest;
        blockIdx+=1;
        size_t last;
        if (lastLast) {
            last=((nFlags>startPos+2)?(nFlags-2):startPos);
        } else {
            last=((nFlags>startPos+1)?(nFlags-1):startPos);
        }
        auto stopBlock=last/bitsPerEl;
        bool wasReset=false;
        while (1){
            while (rest!=0){
                --rest;
                if (atomicCAS(idx,newVal,0)==0){
                    return idx;
                }
                idx+=1;
                val>>=usedBits;
                if (idx==startPos) { // finished search
                    if (opIndex(nFlags-1)==0 && (extraCheck is null || extraCheck[nFlags-1]==0)){
                        if (atomicCAS(nFlags-1,newVal,0)==0){
                            return nFlags-1;
                        }
                    }
                    return nFlags;
                }
                if (idx>last) {
                    if (wasReset) throw new Exception("reset loop",__FILE__,__LINE__);
                    wasReset=true;
                    idx=0;
                    blockIdx=0;
                    stopBlock=startPos/bitsPerEl;
                    if (idx==startPos) { // finished search
                        if (atomicCAS(nFlags-1,newVal,0)==0){
                            return nFlags-1;
                        }
                        return nFlags;
                    }
                    break;
                }
            }
            while(blockIdx!=stopBlock){
                val=data[blockIdx];
                auto r=((((val+add0)|val)^val)&bitToCheck);
                if (r!=0) break;
                blockIdx+=1;
            }
            idx=blockIdx*bitsPerEl;
            rest=bitsPerEl;
            val=(data[blockIdx]);
            blockIdx+=1;
        }
    }
    /// number of flags stored
    size_t length(){
        return nFlags;
    }
    mixin(serializeSome("dchem.FlagsArray","nFlags|data"));
    void desc(CharSink sink){
        auto s=dumper(sink);
        s("{class:\"dchem.FlagsArray\",length:")(nFlags)(",data:");
        auto i=nFlags;
        while(i!=0){
            if(i!=nFlags) s(".");
            --i;
            s(printBits(opIndex(i),2));
        }
        s("}");
    }
    void clearData(){
        data[]=0;
    }
}
