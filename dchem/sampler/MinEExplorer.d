module dchem.sampler.MinEExplorer;
import dchem.Common;
import dchem.sys.ParticleSys;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import tango.math.random.Random;
import blip.rtest.RTest;
import Atomic=blip.sync.Atomic;
import blip.t.math.Math:max;
import tango.util.container.more.Heap;
import blip.container.BatchedGrowableArray;
import tango.util.container.HashSet;
import dchem.input.RootInput;
import blip.container.Deque;

/// unique key for a point, dimensions might change
/// key might have some special values,
///   0: means no point,
///   1: pseudo key for points whose exploration has been excluded (for example due to constraints)
/// the key itself has the first 10 bits (size might change) that are used for collision reduction
/// and the rest that encodes a mostly sequential number (this is indirectly useful as the exploration is not random),
/// and to allow efficient storage of dense arrays (making unbundling of properties of points easy)
/// all keys with 0 at the beginning and only the first 10 bit set are thus never generated,
/// and might be used as special keys (see the special values above)
struct Point{
    ulong data;
    ulong key(){
        return data;
    }
    /// sets base key (should not be too large)
    void key(ulong k){
        assert(k<=0x0000_0FFF_FFFF_FFFFUL,"key too large");
        data=k&0x0000_0FFF_FFFF_FFFFUL; // drop the extra &?
    }
    /// set base key without loosing too many bits (and guranteeing that if it was non 0 it stays non 0)
    static Point opCall(ulong k){
        Point res;
        res.data=(k & 0x0000_0FFF_FFFF_FFFFUL);
        auto r=((k>>21)&0x0000_07FF_FF80_0000UL);
        if (r!=0){
            res.data=(res.data^r)|0x0000_0800_0000_0000UL;
        }
        return res;
    }
    bool isValid(){
        return (data&0x0000_0FFF_FFFF_F000UL)!=0;
    }
    mixin(serializeSome("dchem.Point","data"));
    mixin printOut!();
    hash_t toHash(){
        static if(hash_t.sizeof==4){
            return cast(hash_t)(0xFFFF_FFFF&(data^(data>>32)));
        } else {
            return cast(hash_t)data;
        }
    } 
}

/// unique key for a point, and a direction in the same structure
struct PointAndDir{
    /// 16 bit direction, 36 bit sequence number, 12 bit collision reduction, the bit grouping might change, don't count on it
    ulong data; /// key of a main point (64 bit as trade off between collision probability and being small). 0 is an invalid key
    static PointAndDir opCall(Point k,uint dir){
        PointAndDir res;
        if (dir!=uint.max){
            assert(dir<0xF_FFFF,"direction too large");
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)dir)<<44);
        } else {
            res.data=(k.data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
        return res;
    }
    /// direction
    uint dir(){
        uint res=cast(uint)(data>>44);
        return ((res!=0xF_FFFF)?res:uint.max);
    }
    //sets direction
    void dir(uint d){
        if (d!=uint.max){
            assert(d<0xF_FFFF,"direction too large");
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|((cast(ulong)d)<<44);
        } else {
            data=(data&0x0000_0FFF_FFFF_FFFFUL)|0xFFFF_F000_0000_0000;
        }
    }
    /// base key
    Point point(){
        return Point(0x0000_0FFF_FFFF_FFFFUL & data);
    }
    /// sets base key, dropping extra bits
    void point(Point k){
        data=(0xFFFF_F000_0000_0000 & data)|(k.data & 0x0000_0FFF_FFFF_FFFFUL);
    }
    mixin(serializeSome("dchem.PointAndDir","data"));// split in Point and dir???
    mixin printOut!();
    struct ExpandedPointAndDir{
        PointAndDir pDir;
        uint dir(){
            return pDir.dir();
        }
        void dir(uint d){
            pDir.dir=d;
        }
        Point point(){
            return pDir.point();
        }
        void point(Point p){
            pDir.point=p;
        }
        mixin(serializeSome("dchem.PointDir","point|dir"));// split in Point and dir???
        mixin printOut!();
    }
    ExpandedPointAndDir expanded(){
        ExpandedPointAndDir e;
        e.pDir=*this;
        return e;
    }
}

alias ulong EKey; /// key of a MinEExplorer instance, 0 is an invalid key

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
    size_t findFreeAndSet(size_t startPos,uint newVal=0,bool lastLast=true){
        if (nFlags==0) return 0;
        if (nFlags==1) {
            if (atomicCAS(0,newVal,0)==0){
                return 0;
            }
            return 1;
        }
        assert(startPos<nFlags,"startPos out of bounds");
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
                if (idx==startPos) { // finished search¨
                    if (opIndex(nFlags-1)==0){
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
}

struct LocalPoint(T){
    Point point;
    PointAndDir generator;
    MainPoint!(T) _mainPoint;
    MainPoint!(T) mainPoint(){
        if (_mainPoint !is null) return _mainPoint;
        synchronized{
            if (_mainPoint !is null) return _mainPoint;
            if (! generator.isValid()) throw new Exception(collectAppender(delegate void(CharSink s){
                dumper(s)(point)(" has no mainPoint and generator ")(generator)(" is non valid");
            }),__FILE__,__LINE__);
            
        }
    }
}
/// a main evaluation point
class MainPoint(T){
    enum GFlags{
        None=0,/// no flags
        InProgress=(1<<0),        /// a calculation at this point is in progress
        EvaluatedEnergy=(1<<1),   /// the energy was evaluated
        EvaluatedGradient=(1<<2), /// gradient was evaluated, frame of reference is established
        DoNotExplore=(1<<3),   /// no further exploration should be performed starting from this point
        FullyExplored=(1<<4),  /// all directions around this point have been explored
        FullyEvaluated=(1<<5), /// all directions around this point have been explored and evaluated
        OldApprox=(1<<6),      /// energy/gradient are based on old data, newer better data should be available
        LocalCopy=(1<<7),      /// this MainPoint is a local copy done for efficency reasons (localContext is not the owner)
    }
    MinEExplorer!(T) localContext;
    Point point;
    /// position of the point (energy, derivatives,... might be invalid, only real positions have to be valid)
    ParticleSys!(T) pos;
    /// direction of the minimum in the dual space with norm 1 wrt. euclidische norm
    /// (frame of reference for minimization)
    DynPVect!(T,2) minDir;
    /// directions, at index 0-core Energy 1-core Gradient, index i, direction i/2 of the derivatives in the dual space (DynPVect!(T,2))
    /// (i&1)==0 if positive, 1 if negative, the last index is of direction 0 again.
    /// flags values are the DirFlags: 0-not explored, 1- explored
    FlagsArray dirFlags;
    /// neighbors, and the direction with respect to this point in which they are
    LocalGrowableArray!(PointAndDir) neighbors;
    /// scale of mindir to recover the dual gradient (useful??)
    T minDirScale;
    /// repulsion size (useful to block out regions)
    T repulsionSize;
    /// exploration size for this point (useful to recover positions) 
    T explorationSize;
    /// bit or of GFlags of the current point
    uint gFlags;
    
    /// structure to encode a main point efficienly (useful to transfer it o another context)
    struct DriedPoint{
        Point point;
        ParticleSys!(T) pos;
        DynPVect!(T,2) minDir;
        FlagsArray dirFlagsTopo; /// topological direction flags
        FlagsArray[MKey] dirFlagsMethods;
        LocalGrowableArray!(PointAndDir) neighbors;
        T minDirScale;
        T repulsionSize;
        T explorationSize;
        uint gFlags;
        
    }
    
    /// constructor
    this(MinEExplorer!(T) localContext,Point point,ParticleSys!(T) pos,uint gFlags=GFlags.None,
        T repulsionSize=-1,T explorationSize=-1,DynPVect!(T,2) minDir=null,T minDirScale=T.init,
        FlagsArray dirFlags=null,PointAndDir[] neighbors=null)
    {
        this.localContext=localContext;
        this.point=point;
        this.pos=pos;
        this.gFlags=gFlags;
        this.repulsionSize=repulsionSize;
        this.explorationSize=explorationSize;
        if (this.repulsionSize<0) this.repulsionSize=localContext.repulsionSize;
        if (this.explorationSize<=0) this.explorationSize=localContext.explorationSize;
        this.minDir=minDir;
        this.minDirScale=minDirScale;
        this.dirFlags=dirFlags;
        if (this.dirFlags is null){
            this.dirFlags=new FlagsArray(localContext.nDirs);
        }
        this.neighbors=lGrowableArray(neighbors,neighbors.length,GASharing.Global);
    }
    
    /// explores the next direction, immediately marks it as in exploration
    /// returns an invalid point only if all direction are explored, and an invalid direction but this point if a
    /// frame of exploration is not yet defined
    PointAndDir exploreNext(bool cheapGrad){
        auto exploreFlags=DirFlags.ExploredLocal;
        if ((gFlags&GFlags.DoNotExplore)!=0){
            return PointAndDir(0,uint.max);
        }
        if ((gFlags&GFlags.EvaluatedGradient)==0){
            synchronized(this){
                if ((gFlags&(GFlags.InProgress|GFlags.EvaluatedGradient))==0){
                    gFlags|=GFlags.InProgress;
                    return PointAndDir(point,0);
                }
            }
        }
        if ((gFlags&GFlags.InProgress)!=0){
            return PointAndDir(this.point,uint.max); // needs to wait evaluation
        }
        if ((gFlags&GFlags.FullyExplored)==0){
            synchronized(this){
                if ((gFlags&GFlags.FullyExplored)==0){
                    size_t nextDir=dirFlags.length;
                    if(dirFlags.findFreeAndSet(1,exploreFlags)){
                        nextDir=0;
                    } else {
                        auto nextDir=dirFlags.findFreeAndSet(start,exploreFlags);
                    }
                    if (nextDir!=dirFlags.length) return PointAndDir(this.point,nextDir);
                    gFlags|=GFlags.FullyExplored;
                }
            }
        }
        return PointAndDir(0,0);
    }
    /// returns the dir value for the given dimension and sign
    uint toDir(uint idim,bool neg){
        assert(idim<dirFlags.length/2,"dim out of bounds");
        if (idim==0 && neg) return dirFlags.length-1;
        return 2*idim+1-(neg?1:0);
    }
    /// transforms a non core (i.e. non 0) dir to dimension and sign
    void fromDir(uint dir,out uint dim,out bool neg){
        assert(dir>0,"special core direction not mappable");
        assert(dir<dirFlags.length,"dir is too large");
        dim=dir/2;
        if(dir==dirFlags.length-1) dir=0;
        neg=((dir&1)==0);
    }
    /// direction
    DynPVector!(T,2) dualDir(uint dir){
        if (dir==0) return pos;
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        T val=(neg?-1:1);
        auto v0=minDir.emptyCopy;
        v0[]=0;
        v0[rDir]=val;
        auto newDir=rotateEiV(V,M)(0,minDir,v0);
        return newDir;
    }
    /// calculates the point in the given direction
    ParticleSys!(T) posInDir(uint dir){
        auto newDir=dualDir(dir);
        pos.projectInDualTSpace(newDir);
        auto nAtt=newDir.norm();
        if (nAtt<localContext.minProjectionResidual){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)("discarded because direction is mostly in null space");
            });
            /// try to get rid of all null diretions checking the diagonal of the rotated overlap???
            return null; // continue all the same or return pos???
        }
        /// immediately check expected cartesian norm change and possibly increase length? could be wrong for periodic directions
        newDir*=explorationSize/nAtt;
        
        auto resPSys=pos.dup(PSDupLevel.DynProperties|PSDupLevel.DynPNullDx|PSDupLevel.DynPNullMddx|PSDupLevel.HiddenVars);
        pos.addFromDualTSpace(T)(newDir,resPSys.dynVars.x);
        newDir.giveBack();
        auto err=localContext.constraints.applyR(resPSys);
        if (err>0.1){
            localContext.logMsg(delegate(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)("discarded because we might have a constraint clash");
            });
            return null;
        }
        resPSys.updateHVars();
    }
    /// returns if the given direction is acceptable. If it is accepted sets the local directions covered
    /// and returns a valid Point for it.
    Point acceptNewDirection(DynPVector!(T,0) newPos,uint dir){
        if (dir==0) return point; // check that newPos and pos are very close???
        uint rDir; bool neg;
        fromDir(dir,rDir,neg);
        Point newPoint=Point(0);
        /// verify pos:
        auto diff=newPos.dup();
        diff.axpby(pos.dynVars.x,-1,1);
        if (localContext.minMoveInternal>0 && sqrt(diff.norm2())<localContext.minMoveInternal){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in the internal coordinates of exploration direction ")(dir)(" of point ")
                    (point)(" is too small in internal coordinates ")(sqrt(diff.norm2()))(", discarding evaluation");
            });
            return newPoint;
        }
        auto deriv1=pos.dynVars.emptyDx();
        deriv1[]=0;
        pos.addToTSpace(diff,deriv1);
        diff.giveBack();
        // make an optimized version when the overlap is the identity??
        auto deriv1Dual=pos.toDualTSpace(deriv1,deriv1Dual);
        // original norm in dual space
        auto origNorm=sqrt(deriv1Dual.norm2());
        if (origNorm<localContext.minNormDualSelf()*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is too small:")(origNorm)(", discarding evaluation");
            });
            return newPoint;
        }
        if (origNorm>maxNormDual*explorationSize || origNorm<localContext.minNormDual*explorationSize){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in dual space of exploration in direction ")(dir)(" of point ")(point)
                    (" is outside the exploration zone with length ")(origNorm)("*")(explorationSize)
                    (", continuing evaluation");
            });
        }
        // norm in cartesian units
        normCartesian=sqrt(pos.dotInTSpace(deriv1,deriv1Dual));
        if (normCartesian<localContext.minRealNormSelf){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("norm in real space of exploration in direction ")(dir)(" of point ")(point)
                    (" is too small:")(normCartesian)(", discarding evaluation");
            });
        }
        /// checking which direction of the dual space we are really exploring
        auto deriv2Dual=rotateVEi(minDir,0,deriv1Dual);
        volatile T maxVal=0;
        size_t iMax;
        foreach (i,v;deriv2Dual.pLoop){
            auto vp=abs(v);
            if (vp>maxVal){
                synchronized{
                    if (vp>maxVal){
                        maxVal=vp;
                        iMax=i;
                    }
                }
            }
        }
        auto dirMax=toDir(iMax,deriv2Dual[iMax]<0);
        if (dirMax!=dir){
            if (maxVal>origNorm/sqrt(2)){
                auto actualDirFlags=dirFlags.atomicCAS(iMax,GFlags.Explored,GFlags.Free);
                if (actualDirFlags == GFlags.Free){
                    localContext.logMsg(delegate void(CharSink s){
                        dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                            (" is not mostly in the expected direction in the dual space, but in dir ")
                            (iMax)(", continuing evaluation declaring also it as visited (check distances also???)");
                    });
                } else {
                    localContext.logMsg(delegate void(CharSink s){
                        dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                            (" is not mostly in the expected direction in the dual space, but in dir ")
                            (iMax)(", and that was already visited, discarding evaluation");
                    });
                    return newPoint;
                }
            } else {
                localContext.logMsg(delegate void(CharSink s){
                    dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                        (" is not mostly in the expected direction in the dual space, but in dir ")
                        (iMax)(", but with angle larger than 45°, continuing evaluation");
                });
            }
        }        
        auto deriv2=rotateVEi(minDir,0,deriv1);
        deriv2Dual.giveBack();
        // length of the directions in cartesian units
        auto ov=pos.maybeOverlap();
        scope rotOv=rotateVEi(minDir,0,ov.data.T.dup(true));
        auto lenDirs=diag(rotOv);
        unaryOpStr!("*aPtr0=sqrt(*aPtr0);")(lenDirs);
        
        // adds new point
        newPoint=localContext.nextPointId();
        localContext.addLocalPoint(new MainPoint!(T)(localContext,newPoint,newPos));
        
        PointAndDir[128] buf;
        auto neighAtt=lGrowableArray(buf,0);
        neighAtt.append(PointAndDir(newPoint,dir));
        if (dirMax!=dir)
            neighAtt.append(PointAndDir(newPoint,dirMax));
        
        // check that we are still mostly in the dir direction in the cartesian metric
        if (abs(deriv2[rdir])<sameDirCosAngle*normCartesian*lenDirs[rdir] || (deriv2[rdir]<0 && (dir&1)==1)){
            localContext.logMsg(delegate void(CharSink s){
                dumper(s)("exploration in direction ")(dir)(" of point ")(point)
                    (" is not mostly in the expected direction in the cartesian space: cos(angle):")
                    (abs(deriv2[rdir])/normCartesian*lenDirs[rdir])(", continuing evaluation");
            });
        }
        // check other directions that we might cover in cartesian units
        auto ndir=dirFlags.length/2;
        for(idir=0;idir<ndir;++idir){
            if (abs(deriv2[idir])>=sameDirCosAngle*explorationSize*lenDirs[i]){
                auto dirAtt=toDir(idir,deriv2[idir]<0);
                if (dirAtt != dir && dirAtt!= iMaxDir){
                    dirFlags.atomicCAS(dirAtt,DirFlags.Explored,DirFlags.Free);
                    if (deriv2[idir]!=0){ // don't store directions for invalid dirs
                        neighAtt.append(PointAndDir(newPoint,dirAtt));
                    }
                }
            }
        }
        synchronized(this){
            neighbors.appendArr(neighAtt.data); // duplicates are stored contiguously, first el is original search dir
        }
        return newPoint;
    }
    /// the various kinds of neighbors
    enum NeighKind{
        NoNeigh, /// not a neighbor
        NeighInDir, /// a neighbor that explores one or more directions of this point
        NeighNoDir, /// a neighbor that does not explore a direction of this point
        NeighVeryClose, /// a point very close to this point
    }
    /// checks if the point passed is a neighbor of this point, if it is checks the directions blocked by this
    /// point. Returns the kind of neighbor that the new point is
    NeighKind checkNeighbors(DynPVector!(T,0)newPos,Point p) {
        assert(0,"to do");
        return NeighKind.NoNeigh;
    }
    bool evalWithContext(CalculationContext c){
        bool calcE=(gFlags&GFlags.EvaluatedEnergy)==0;
        bool calcF=((gFlags&GFlags.EvaluatedGradient)==0 || cheapGrad());
        if (calcE||calcF){
            Real e=pSys.dynVars.potentialEnergy;
            mixin(withPSys(`pSys[]=pos;`,"c."));
            c.updateEF(calcE,calcF);
            mixin(withPSys(`pos[]=pSys;`,"c."));
            if (isNAN(pos.dynVars.potentialEnergy)) {
                pos.dynVars.potentialEnergy=e; // replace with if(!calcE) ?
            } else {
                e=pos.dynVars.potentialEnergy;
            }
            if (!calcF){
                if (pos.potentialEnergy<=localContext.toExploreMore.peek.energy){
                    calcF=true;
                    c.updateEF(false,true);
                    if (isNAN(pos.dynVars.potentialEnergy)) {
                        pos.dynVars.potentialEnergy=e; // replace with if(!calcE) ?
                    } else {
                        e=pos.dynVars.potentialEnergy;
                    }
                }
            }
            auto target=localContext.activeExplorers[localContext.owner];
            if (calcF){
                if ((gFlags&GFlags.LocalCopy)!=0){
                    assert(target is activeExplorers[localContext.key],"non local copy and not local");
                    target.addGradEvalLocal(point,pSysWriter(pos));
                } else {
                    synchronized(this){
                        if (calcF) {
                            localContext.didGradEval(localContext.key,point);
                            gFlags|=GFlags.EvaluatedEnergy|GFlags.EvaluatedGradient;
                        } else {
                            gFlags|=GFlags.EvaluatedEnergy;
                        }
                    }
                }
            } else { // just updated energy
                target.addEnergyEvalLocal(point,e);
            }
            return true;
        } else {
            return false;
        }
    }
}

// object that keeps the journal of the computations done
class MinEJournal{
    char[] journalName;
    BinSink jSink;
    SBinSerializer jSerial;
    const char[] jVersion="MinEJournal v1.0";
    struct MainPointAdd(T){
        Point point;
        PSysWriter!(T) pos;
        T stepSize;
        T repulsionSize;
        
        void reapply(MinEExplorer expl){}
    }
    struct EFEval(T){
        Point key;
        ParticleSys!(T) pos;
        Real energy;

        void reapply(MinEExplorer expl){}
    }
    struct EEval(T){
        Point!(T) p;
        Real energy;

        void reapply(MinEExplorer expl){}
    }
    struct NeighAdd(T){
        Point!(T) p1;
        Point!(T) p2;

        void reapply(MinEExplorer expl){}
    }
    
}

class MinHeapSync(T){
    MinHeap!(T) heap;
    WaitCondition nonEmpty;
    
    bool nonEmptyHeap(){
        return heap.length!=0;
    }
    
    this(){
        nonEmpty=new WaitCondition(&nonEmptyHeap);
    }
    void push(T[] t){
        synchronized(this){
            heap.push(t);
        }
        nonEmpty.checkCondition();
    }
    void push(T t){
        synchronized(this){
            heap.push(t);
        }
        nonEmpty.checkCondition();
    }
    T pop(){
        synchronized(this){
            return heap.pop();
        }
    }
    /// returns the minimal energy elements, waits if no elements is available until some becomese available
    T waitPop(){
        while (1){
            synchronized(this){
                if (heap.length>0)
                    return heap.pop();
            }
            nonEmpty.wait();
        }
    }
}

interface ActiveExplorer(T){
    /// stops all explorers (at the moment there is no support for dynamic adding/removal of explorers, worth adding???)
    void shutdown(EKey);
    /// key of this explorer
    EKey key();
    /// where the journal is kept
    char[] journalPos();
    /// load in some arbitrary units
    Real load();
    
    /// adds energy (and removes from the inProgress) for a local point, then bCasts addEnergyEval
    void addEnergyEvalLocal(Point p,Real energy);
    /// adds energy to the evaluated ones
    void addEnergyEval(Point p,Real energy);
    /// adds gradient value to a point that should be owned. Energy if not NAN replaces the previous value
    void addGradEvalLocal(Point p,PSysWriter!(T) pSys);
    /// communicates that the given point is being expored
    /// flags: communicate doubleEval?
    void addExploredPoint(EKey owner,Point point,PSysWriter!(T) pos,uint flags);
    
    void finishedExploringPointLocal(Point);
    /// finished exploring the given point (remove it from the active points)
    void finishedExploringPoint(Point);
    /// drops all calculation/storage connected with the given point, the point will be added with another key
    /// (called upon collisions)
    void dropPoint(EKey,Point);
    
    // /// if the given point is still active
    // bool pointIsExplorable(EKey,Point);
    
    /// returns the position of the given point
    PSysWriter!(T)pointPos(Point);
    /// energy for the point (NAN if not yet known)
    Real energyForPoint(Point);
    /// returns a point to evaluate
    Point pointToEvaluate();
    /// called when an evaluation fails, flags: attemptRetry/don't Retry
    void evaluationFailed(Point,uint flags);
}
alias HashSet Set;
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

/// a structure to keep point and its energy together
struct PointAndEnergy{
    Point point;
    Real energy;
    mixin(serializeSome("PointAndEnergy","point|energy"));
    mixin printOut!();
}

struct LoadStats{
    Real load;
    EKey explorer;
    mixin(serializeSome("LoadStats","load|explorer"));
    mixin printOut!();
}

class RemotePointEval(T):RemoteTask{
    Point point;
    char[] ownerUrl;
    CalculationContext ctx;
    mixin(serializeSome("dchem.minEE.RemotePointEval!("~T.stringof~")","point|ownerUrl"));
    
    this(){}
    this(Point p,char[]oUrl){
        point=p;
        ownerUrl=oUrl;
    }
    void execute(Variant args){
        auto ctx=cast(ActiveExplorer!(T))cast(Object)ProtocolHandler.proxyForUrl(ownerUrl);
        auto pPos=ctx.pointPos(point); // get real point
        ctx=args.get!(CalculationContext)();
    }
    void stop(){
        if (ctx!is null) ctx.stop();
        ctx=null;
    }
    
}
class MasterCalculator{
    Set!(PointAndDir) inExploration; /// point directions that are in exploration
    Set!(Point) inEvaluation; /// points that are in evaluation
    Set!(Point) toEvaluate; /// points that should be avaluated
    Deque!(CalculationContext) waitingContexts;
    size_t overPrepare=0;
    /// adds a context that can calculate something
    void addContext(CalculationContext c){
        synchronized(this){
            waitingContexts.append(c);
        }
        update();
    }
    /// adds a point to evaluate coming from PointAndDir
    void addPointForPointDir(Point p,PointAndDir origin){
        synchronized(this){
            if (p.isValid){
                toEvaluate.add(p);
            }
            inExploration.remove(origin);
        }
        update();
    }
    /// updates the calculator: tries to give work to waitingContexts and
    void update(){
        Point p;
        CalculationContext ctx;
        while(1){
            ctx=null;
            synchronized(this){
                if (!toEvaluate.take(p)) break;
                if (!waitingContexts.popFront(ctx)) break;
            }
            // ctx.remoteExe...
        }
        if (p.isValid && ctx is null) {
            synchronized(this){
                toEvaluate.add(p);
            }
        }
    }
}

class MinEExplorer(T): Sampler,ActiveExplorer!(T){
    static UniqueNumber!(ulong) nextLocalId;
    static UniqueNumber!(ulong) nextPntNr;
    char[] nameId;
    EKey key;
    /// points to explore more ordered by energy (at the moment this is replicated)
    /// this could be either distribued, or limited to a given max size + refill upon request
    MinHeapSync!(PointAndEnergy) toExploreMore;
    /// list of all the currently active explorers (including this one)
    ActiveExplorer!(T)[EKey] activeExplorers;
    MinHeapSync!(LoadStats) loads; 
    BatchGrowableArray!(Point,batchSize) localPointsKeys; // indexes for local owned points
    HashMap!(Point,EKey) owner;
    CalculationInstance[Point] localCalcInProgress;
    Set!(Point) calcInProgress;
    DenseLocalPointArray!(MainPoint!(T)) localPoints; /// local points (position,...)
    Random randLocal;
    RandomSync rand;
    /// step used for the discretization
    T discretizationStep;
    /// place to store the trajectory (journal)
    char[] trajDir;
    /// the current method for energy evaluation
    InputField evaluator;
    /// constraints (taken from the evaluator)
    MultiConstraint _constraints;
    /// place to log messages
    CharSink log;
    /// notification center
    NotificationCenter nCenter;
    /// unit of measurement for 
    
    /// adds the computed gradient to the point p
    void addGradEval(Point p,PSysWriter!(T) pSys){
        auto mp=localPoints[p];
        auto e=mp.pos.dynVars.potentialEnergy;
        synchronized(mp){
            mp.pos[]=pSys;
            if (isNAN(mp.pos.dynVars.potentialEnergy)) mp.pos.dynVars.potentialEnergy=e;
            gFlags|=GFlags.EvaluatedEnergy|GFlags.EvaluatedGradient;
        }
    }
    /// writes out a log message
    void logMsg(void delegate(CharSink)writer){
        char[512] buf;
        auto gArr=lGrowableArray(buf,0);
        auto s=dumper(&gArr.appendArr);
        s("<MinEELog id=\"")(key)("\" time=\"")(toString(Clock.now()))("\">");
        writer(s.call);
        s("</MinEELog>\n");
    }
    /// writes out a log message
    void logMsg(char[]msg){
        logMsg(delegate void(CharSink s){ s(msg); });
    }
    /// constraints of the current system
    MultiConstraint constraints(){
        return _constraints;
    }
    /// if the gradient is cheap to compute
    bool cheapGrad(){
        return false;
    }
    /// minimum norm of uint vector after projection of null direction to be considered as a valid direction
    T minProjectionResidual(){
        return 0.01;
    }
    /// radius in which to look for neighbors (in units of explorationSize)
    T sameDirCosAngle(){
        return 0.866;
    }
    /// minimum norm in the dual T space to accept an exploration (in units of explorationSize)
    T minNormDual(){
        return 0.4;
    }
    /// minimum norm in the dual T space to accept an exploration for a self generated direction 
    /// (in units of explorationSize)
    T minNormDualSelf(){
        return 0.1;
    }
    /// minimum norm in the real (cartesian) space for a self generated direction to be accept
    T minRealNormSelf(){
        return 0.05;
    }
    /// max norm in the dual T space to accept an exploration (in units of explorationSize)
    T maxNormDual(){
        return 1.9;
    }
    /// explorationSize
    T explorationSize(){ return discretizationStep; }
    /// repulsionSize
    T repulsionSize(){ return discretizationStep; }
    bool pointIsExplorable(EKey eK,Point p){
        MainPoint!(T) mainP;
        synchronized(this){ // sync could be avoided most of the time 
            mainP=localPoints[p];
        }
        if (mainP is null) throw new Exception("asked unknown point "~p.toString(),__FILE__,__LINE__);
        auto flags=mainP.gFlags;
        return (flags&GFlags.Evaluated)&&(!(flags&(GFlags.DoNotExplore|GFlags.FullyExplored|GFlags.FullyEvaluated|GFlags.OldApprox)));
    }
    /// returns a globally unique string 
    char[] nextUniqueStr(){
        return collectAppender(delegate void(CharSinker s){
            s(nameId); s("_"); writeOut(s,nextUniqueId.next());
        });
    }
    /// returns a most likely valid point id
    Point nextPointId(){
        ushort r;
        rand(r);
        return Point(((nextPntNr.next())<<12)|cast(ulong)(r&0xFFF));
    }
    
    /// returns true if the evaluation of the given point is in progress
    bool isInProgress(Point p){
        bool res=false;
        synchronized(this){
            res=(p in calcInProgress)!is null;
        }
        return res;
    }
    
    mixin(serializeSome("dchem.MinEExplorer_"~T.stringof,
        `trajDir: directory where to store the trajectory (journal)`));
    mixin printOut!();
    
    void run(){
        bool restarted=false;
        evaluator.method.setupCalculatorClass();
        // possiby restarts
        if (! restarted){
            cInstance=getInstanceForClass(InstanceGetFlags.ReuseCache|InstanceGetFlags.NoAllocSubOpt|InstanceGetFlags.Wait);
        }
        master.nextCalculation();
    }
    
    // exploration:
    // find smallest energy still "free"
    // evaluate next direction of it
    // apply constraints, if movement is too small declare it as fully visited and neighbor 1
    // bcast point as explored
    // possibly wait for non collision confirmation
    // start calculation
    // if too close to existing points stop calculation???
    // when calculation is finished
    // if mainpoint:
    //    gradient -> orient neighbors, compile visited flags, perform topology analysis (attractor, min,max,...)
    //    second deriv check
    // else store energy, first deriv check??
    
    /// evaluates the next computation
    PointAndDir nextComputation(){
        auto smallPoint=toExploreMore.waitPop();
        auto nextPDir=activeExplorers[owner[smallPoint.point]].exploreNext(cheapGrad);
        if (nextPDir.isValid && nextPDir.dir!=0){
            toExploreMore.push(smallPoint);
        } // else other explorers have a point that is not executable as first, but that should create no problems
        return nextPDir;
    }
    
    /// returns the key of the worker with the lowest expected load
    EKey findNextWorker(){
        auto w=loads.popWait();
        mixin(mkActionMixin("updLoad","w|activeExplorers|loads",`
        w.load=activeExplorers[w.key].load();
        loads.push(w);`));
        if (key==w.key){ // local update
            updLoad();
        } else {
            Task("updLoad",&updLoad).autorelease.submitYield();
        }
        return w.key;
    }
    /// evaluates the given point and dir, and either request a real evaluation there, or communicates that it is blocked
    Point evaluatePointAndDir(){
        activeExplorers[findNextWorker()];
    }
    
    void stop(){
        
    }
    
    bool verify(CharSink s){ return true; }
}
