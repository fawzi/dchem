module dchem.neigh.HierarchicalSort;
import dchem.Common;
import dchem.neigh.NeighModels;
import blip.util.RefCount;
import blip.parallel.smp.PLoopHelpers;

/// 10 bits x 3 index, 2 bits unused
enum TreeMask:uint{
    Idx=0x3FFF_FFFF,
    Pos=0x3FFF_FFFF,
    Rest=~Idx;
    IdxX=0x3FF0_0000,
    IdxY=0x000F_FC00,
    IdxZ=0x0000_03FF,
    IdxBase=0x3FF,
    Unused=0xC000_0000,
}

enum TreeShifts:int{
    Idx=0,
    IdxX=20,
    IdxY=10,
    IdxZ=0,
}

enum TreeConsts{
    MaxSub=10,
}

struct TreePos{
    uint data;
    
    equals_t opEquals(TreeIdx o){
        return data==o.data;
    }
    int opCmp(TreeIdx o){
        return ((data<o.data)?-1:((data==o.data)?0:1));
    }
    
    static TreePos opCall(TreeIdx p,int subDiv=TreeConsts.MaxSub){
        uint idxX=cast(uint)(TreeMask.IdxBase&(p.data>>TreeShifts.IdxX));
        uint idxY=cast(uint)(TreeMask.IdxBase&(p.data>>TreeShifts.IdxY));
        uint idxZ=cast(uint)(TreeMask.IdxBase&(p.data>>TreeShifts.IdxZ));
        return opCall(idxX,idxY,idxZ,subDiv);
    }
    static TreePos opCall(uint idxX, uint idxY, uint idxZ, int subDiv=TreeConsts.MaxSub){
        assert(idxX<=TreeMask.IdxBase,"idxX out of range");
        assert(idxY<=TreeMask.IdxBase,"idxY out of range");
        assert(idxZ<=TreeMask.IdxBase,"idxZ out of range");
        assert(subDiv>0 && subDiv<=TreeConsts.MaxSub,"subdiv out of range");
        TreePos res;
        const uint[] pX=[0b000_000_000,0b000_000_100,0b000_100_000,
            0b000_100_100,0b100_000_000,0b100_000_100,
            0b100_100_100];
        const uint[] pY=[0b000_000_000,0b000_000_010,0b000_010_000,
            0b000_010_010,0b010_000_000,0b010_000_010,
            0b010_010_010];
        const uint[] pZ=[0b000_000_000,0b000_000_001,0b000_001_000,
            0b000_001_001,0b001_000_000,0b001_000_001,
            0b001_001_001];
        int shiftNow=TreeConsts.MaxSub-subDiv;
        idxX>>=shiftNow;
        idxY>>=shiftNow;
        idxZ>>=shiftNow;
        while(shiftNow<TreeConsts.MaxSub){
            pos+=(pX[(idxX&7)]+pY[(idxY&7)]+pZ[(idxZ&7)])<<shiftNow;
            idxX>>=3;
            idxY>>=3;
            idxZ>>=3;
            shiftNow+=3;
        }
        res.data=pos;
        return res;
    }
}

struct TreeIdxT(T){
    
    /// 10 bits x 3 index, 2 bits unused
    ulong data;
    
    equals_t opEquals(TreeIdx o){
        return data==o.data;
    }
    int opCmp(TreeIdx o){
        return ((data<o.data)?-1:((data==o.data)?0:1));
    }
    
    static TreeIdx opCall(uint idxX,uint idxY,uint idxZ){
        TreeIdx res;
        res.data=(idxX<<TreeShifts.IdxX)
            +(idxY<<TreeShifts.IdxY)
            +(idxZ<<TreeShifts.IdxZ);
        return res;
    }
    
    static TreeIdx opCall(TreePos p,int subDiv=TreeConsts.MaxSub){
        assert(subDiv>0&&subDiv<=TreeConsts.MaxSub,"subDiv out of range");
        auto d=(p.data>>TreeShifts.Idx);
        uint idxX=0,idxY=0,idxZ=0;
        int shiftNow=TreeConsts.MaxSub-subDiv;
        d>>=3*shiftNow;
        while(shiftNow<TreeConsts.MaxSub){
            switch(d&0b100_100_100){
            case 0b000_000_000:
                idxX+=0b000<<shiftNow;
                break;
            case 0b000_000_100:
                idxX+=0b001<<shiftNow;
                break;
            case 0b000_100_000:
                idxX+=0b010<<shiftNow;
                break;
            case 0b000_100_100:
                idxX+=0b011<<shiftNow;
                break;
            case 0b100_000_000:
                idxX+=0b100<<shiftNow;
                break;
            case 0b100_000_100:
                idxX+=0b101<<shiftNow;
                break;
            case 0b100_100_000:
                idxX+=0b110<<shiftNow;
                break;
            case 0b100_100_100:
                idxX+=0b111<<shiftNow;
                break;
            default:
                assert(0);
            }
            switch(d&0b010_010_010){
            case 0b000_000_000:
                idxY+=0b000<<shiftNow;
                break;
            case 0b000_000_010:
                idxY+=0b001<<shiftNow;
                break;
            case 0b000_010_000:
                idxY+=0b010<<shiftNow;
                break;
            case 0b000_010_010:
                idxY+=0b011<<shiftNow;
                break;
            case 0b010_000_000:
                idxY+=0b100<<shiftNow;
                break;
            case 0b010_000_010:
                idxY+=0b101<<shiftNow;
                break;
            case 0b010_010_000:
                idxY+=0b110<<shiftNow;
                break;
            case 0b010_010_010:
                idxY+=0b111<<shiftNow;
                break;
            default:
                assert(0);
            }
            switch(d&0b001_001_001){
            case 0b000_000_000:
                idxZ+=0b000<<shiftNow;
                break;
            case 0b000_000_001:
                idxZ+=0b001<<shiftNow;
                break;
            case 0b000_001_000:
                idxZ+=0b010<<shiftNow;
                break;
            case 0b000_001_001:
                idxZ+=0b011<<shiftNow;
                break;
            case 0b001_000_000:
                idxZ+=0b100<<shiftNow;
                break;
            case 0b001_000_001:
                idxZ+=0b101<<shiftNow;
                break;
            case 0b001_001_000:
                idxZ+=0b110<<shiftNow;
                break;
            case 0b001_001_001:
                idxZ+=0b111<<shiftNow;
                break;
            default:
                assert(0);
            }
            d>>=9;
            shiftNow+=3;
        }
    }
}

struct HPoint(T){
    TreePos treePos;
    TreeIdx treeIdx;
    PIndex particle;
    index_type iPos;
    CellShift cellShift;
    Vector!(T,3) pos;
    equals_t opEquals(TreeIdx o){
        return treeIdx==o.treeIdx && particle==o.particle && iPos==o.iPos && cellShift==o.cellShift;
    }
    int opCmp(ref TreeIdx o){
        if (&o==this) return 0;
        auto k=particle.kind();
        auto k2=o.particle.kind();
        if (k<k2) return -1;
        if (k>k2) return 1;
        int c=treePos.opCmp(o.treePos);
        if (c==0) {
            c=particle.opCmp(o.particle);
            if (c==0){
                c=((iPos<o.iPos)?-1:((iPos==o.iPos)?0:1));
            }
            if (c==0){
                c=((cellShift.data<o.cellShift.data)?-1:((cellShift.data==o.cellShift.data)?0:1));
            }
        }
        return c;
    }
    static HPoint opCall(TreeIdx treeIdx,TreePos treePos,PIndex particle,index_type iPos,CellShift cellShift,Vector!(T,3)pos){
        HPoint res;
        res.treeIdx=treeIdx;
        res.treePos=treePos;
        res.particle=particle;
        res.iPos=iPos;
        res.cellShift=cellShift;
        res.pos=pos;
        return res;
    }
}

class HierarchicalSort(T):NeighList!(T){
    alias T dtype;
    alias SafeType!(Vector!(T,3)) SafeV;
    
    Cell!(T) cell;
    SegmentedArray!(Vector!(T,3)) pos;
    T scale,minDist,maxDist;
    Vector!(T,3) center;
    BulkArray!(TreePos) sortedTreePos; /// use this more compact array for binary search
    BulkArray!(HPoint!(T)) sortedHPoints; /// all info on the points
    SegmentedArray!(size_t) indexOfPos; /// backref particle,index->sorted el
    size_t[] kindStarts;
    KindRange  kRange;
    int levels;
    
    this(Cell!(T) cell,SegmentedArray!(Vector!(T,3)) pos,T minGridsize=0.4){
        this.kRange=pos.kRange;
        this.cell=cell;
        this.pos=pos;
        Vector!(T,3) cm;
        foreach(i,pk,lk,r;pos.sLoop){
            cm+=r;
        }
        cm*=(cast(T)1)/(cast(T)pos.length);
        auto period=cell.periodicFlags;
        Vector!(T,3)[3] periodDirs;
        Vector!(T,3)[3] periodDirsN;
        int nPeriod=0;
        if ((period&CellPeriodic.x)!=0){
            periodDirs[0]=Vector!(T,3)(cell.h[0,0],cell.h[1,0],cell.h[2,0]);
            nPeriod=1;
        }
        if ((period&CellPeriodic.y)!=0){
            periodDirs[nPeriod]=Vector!(T,3)(cell.h[0,1],cell.h[1,1],cell.h[2,1]);
            ++nPeriod;
        }
        if ((period&CellPeriodic.z)!=0){
            periodDirs[nPeriod]=Vector!(T,3)(cell.h[0,2],cell.h[1,2],cell.h[2,2]);
            ++nPeriod;
        }
        for (int i=1;i<nPeriod;++i){
            periodDirsN[i]=periodDirs[i]*((cast(T)1)/periodDirs.norm2);
            for (int j=0;j<i;++j){
                periodDirsN[i]-=dot(periodDirsN[i],periodDirsN[j])*periodDirsN[j];
            }
            cm-=dot(cm,periodDirsN[i])*periodDirsN[i];
        }
        minDist=0;
        if (period!=CellPeriodic.xyz){
            T dMax2=0;
            foreach(i,pk,lk,r;pos.sLoop){
                auto rAtt=r-cm;
                for(int j=0;j<nPeriod;++j){
                    rAtt-=dot(periodDirsN[j],rAtt)*periodDirsN[j];
                }
                auto dAtt2=rAtt.norm22;
                if (dMax2<dAtt2) dMax2=dAtt2;
            }
            scale=sqrt(dMax2);
            cm+=scale;
        } else {
            scale=0;
        }
        if (period!=CellPeriodic.None){
            scope nn=zeros([nPeriod,nPeriod]);
            T dMax2=0;
            for (int i=0;i<nPeriod;++i)
            for (int j=0;j<i;++j){
                auto dNow2=dot(periodDirs[i],periodDirs[j]);
                nn[i,j]=dNow2;
                nn[j,i]=dNow2;
                if (dMax2<dNow) dMax2=dNow2;
            }
            NArray!(ComplexTypeOf!(T),1) eigV=eig(nn);
            T minEv=eigV[0].re,maxEv=eigV[0].re;
            for (int i=1;i<3;++i){
                if (eigV[i].re<minEv) minEv=eigV[i].re;
                if (eigV[i].re>maxEv) maxEv=eigV[i].re;
            }
            minDist=0.5*sqrt(abs(minEv));
            maxDist=0.5*sqrt(abs(maxEv));
            scale=max(scale,3*dMax2);
        } else {
            minDist=T.max;
            maxDist=T.max;
        }
        scale+=cast(T)1/cast(T)1000;
        this.levels=scale/minGridsize;
        this.center=cm;
        auto blowupFactor=1;
        if ((period&CellPeriodic.x)!=0) blowupFactor*=2;
        if ((period&CellPeriodic.y)!=0) blowupFactor*=2;
        if ((period&CellPeriodic.z)!=0) blowupFactor*=2;
        this.sortedHPoints=BulkArray!(HPoint!(T))(blowupFactor*pos.length);
        this.sortedTreePos=BulkArray!(TreePos)(blowupFactor*pos.length);
        auto aStruct=pos.arrayStruct;
        this.kindStarts=new size_t[](kRange.length+1);
        KindIdx sStart=aStruct.kRange.kStart;
        auto refPos=aStruct.kindStarts[kRange.kStart-sStart];
        foreach(i,k;kRange){
            this.sortedTreePos[i]=(aStruct.kindStarts[k-sStart]-refPos)*blowupFactor;
        }
        size_t iPoint=0;
        foreach(iPos,pk,lk,v,r;pos.sLoop){
            auto redPos=cell.h_inv*r;
            if ((period&CellPeriodic.x)!=0) redPos.x=redPos.x-floor(redPos.x);
            if ((period&CellPeriodic.y)!=0) redPos.y=redPos.y-floor(redPos.y);
            if ((period&CellPeriodic.z)!=0) redPos.z=redPos.z-floor(redPos.z);
            redPos-=this.center;
            Vector!(SafeType!(T),3)[8] posAtt;
            CellShift[8] cellShift;
            int nPos=1;
            posAtt[0]=cellShift.h*p;
            cellShift[0]=CellShift(0,0,0);
            if ((period&CellPeriodic.x)!=0){
                for(int iP=0;iP<nPos;++iP){
                    cellShift[nPos+iP]=cellShift[iP];
                    if (redPos[0]<cast(T)1/cast(T)2){
                        posAtt[nPos+iP]=posAtt[iP]+Vector!(T,3)(cell.h[0,0],cell.h[1,0],cell.h[2,0]);
                        cellShift[nPos+iP].xShift=1;
                    } else {
                        posAtt[nPos+iP]=posAtt[iP]-Vector!(T,3)(cell.h[0,0],cell.h[1,0],cell.h[2,0]);
                        cellShift[nPos+iP].xShift=-1;
                    }
                nPos*=2;
            }
            if ((period&CellPeriodic.y)!=0){
                for(int iP=0;iP<nPos;++iP){
                    cellShift[nPos+iP]=cellShift[iP];
                    if (redPos[1]<cast(T)1/cast(T)2){
                        posAtt[nPos+iP]=posAtt[iP]+Vector!(T,3)(cell.h[0,1],cell.h[1,1],cell.h[2,1]);
                        cellShift[nPos+iP].yShift=1;
                    } else {
                        posAtt[nPos+iP]=posAtt[iP]-Vector!(T,3)(cell.h[0,1],cell.h[1,1],cell.h[2,1]);
                        cellShift[nPos+iP].yShift=-1;
                    }
                }
                nPos*=2;
            }
            if ((period&CellPeriodic.z)!=0){
                for(int iP=0;iP<nPos;++iP){
                    cellShift[nPos+iP]=cellShift[iP];
                    if (redPos[0]<cast(T)1/cast(T)2){
                        posAtt[nPos+iP]=posAtt[iP]+Vector!(T,3)(cell.h[0,2],cell.h[1,2],cell.h[2,2]);
                        cellShift[nPos+iP].zShift=1;
                    } else {
                        posAtt[nPos+iP]=posAtt[iP]-Vector!(T,3)(cell.h[0,2],cell.h[1,2],cell.h[2,2]);
                        cellShift[nPos+iP].zShift=-1;
                    }
                }
                nPos*=2;
            }
            for (int iP=0;iP<nPos;++iP){
                HPoint!(T)*hAtt=sortedHPoints.ptrI(iPoint);
                ++iPoint;
                hAtt.pos=posAtt[iP];
                hAtt.cellShift=cellShift[iP];
                hAtt.particle=pk;
                hAtt.iPos=iPos;
                TreeIdx treeIdx;
                
                auto pX=posAtt[iP].x/scale+(cast(T)1/cast(T)2);
                assert(pX>0&&pX<1,"position out of scale");
                auto iValX=cast(uint)(pX*uint.max);
                iValX>>=(cast(int)typeof(iValX).sizeof*4-TreeConsts.MaxSub);
                auto pY=posAtt[iP].y/scale+(cast(T)1/cast(T)2);
                assert(pY>0&&pY<1,"position out of scale");
                auto iValY=cast(uint)(pY*uint.max);
                iValY>>=(cast(int)typeof(iValY).sizeof*4-TreeConsts.MaxSub);
                auto pZ=posAtt[iP].z/scale+(cast(T)1/cast(T)2);
                assert(pZ>0&&pZ<1,"position out of scale");
                auto iValZ=cast(uint)(pZ*uint.max);
                iValZ>>=(cast(int)typeof(iValZ).sizeof*4-TreeConsts.MaxSub);
                // mask out non used bits (levels..TreeConsts.MaxSub+1) in iVal*??
                auto treeIdx=TreeIdx(iValX,iValY,iValZ);
                hAtt.treeIdx=treeIdx;
                hAtt.treePos=TreePos(treeIdx,levels);
            }
        }
        sortedHPoints.sort(); // compare could be simplified if we would sort kind by kind
        indexOfPos=new SegmentedArray!(size_t)(new SegArrMemMap!(size_t)(pos.arrayStruct));
        foreach (i,ref hAtt;sortedPoints){
            sortedTreePos[i]=hAtt.treePos;
            if(hAtt.cellShift.data==0){
                indexOfPos[hAtt.particle,hAtt.iPos]=i;
            }
        }
    }
    /// a radius, both in real value and its integer value, and magnitude (shift)
    static struct RVal{
        T r; /// radius
        SafeT cutoff2; /// safe radius squared (used as cutoff)
        int shift=-1; /// magnitude (highest bit of iR, but at least the smallest discrete radius)
        uint iR; /// integer radius
        bool valid(){
            return shift>=0;
        }
        static RVal opCall(T r,int levels){
            RVal res;
            res.r=r;
            res.cutoff2=pow2(toSafeType(r));
            if (r>minDist){
                return res;
            }
            auto rVal=(r/scala)*uint.max;
            auto iRval=cast(uint)(rVal*uint.max);
            iRval>>=(cast(int)typeof(iRval).sizeof*4-TreeConsts.MaxSub);
            int shift;
            if (iRval==0){
                shift=-1;
            } else {
                shift=bsr(iRval);
            }
            int minShift=TreeConst.MaxLevels-levels;
            if (shift<minShift){
                shift=minShift;
            } else if (shift>=TreeConst.MaxLevels) {
                shift=-1;
            }
            return res;
        }
    }
    
    /// finds the TreePos ranges that are within r of posAtt
    TreePos[2][] rangesR(T)(TreeIdx posAtt,RVal r,TreePos[2][] ranges){
        assert(r.valid(),"rangesR called with non valid r");
        uint iRval=r.iR;
        int shiftR=r.shift;
        uint lowBit=(1u<<(shiftR+1);
        uint maskR=(1u<<(shiftR+1))-1;
        uint maskR2=(1u<<(shiftR+2))-1;
        uint notMaskR= ~maskR;
        uint[3] iVals;
        iVals[0]=posAtt.x;
        iVals[1]=posAtt.y;
        iVals[2]=posAtt.z;
        uint[4][3] tVals;
        int[3] nIVals=[1,1,1];
        for (int idir=0;idir<3;++idir){
            auto v=iVals[idir];
            uint  low,high;
            if (v<iRval){
                low=0;
            } else {
                low=((v-iRval)&(~maskR));
            }
            if (v+iRval>TreeConst.MaxVal){
                high=TreeConst.MaxVal;
            } else {
                high=((v+iRval)|maskR);
            }
            if (((low^high)&maskR2)!=0){
                // split
                tVals[idir][0]=low;
                tVals[idir][1]=(low|maskR);
                tVals[idir][2]=(high&(~maskR));
                tVals[idir][3]=high;
                nIVals[idir]=2;
            } else {
                tVals[idir][0]=low;
                tVals[idir][1]=high;
                nIVals[idir]=1;
            }
        }
        if (ranges.length<8) ranges=new ranges[2][8];
        int ii=0;
        for(ix=0;ix<nIVals[0];++ix)
        for(iy=0;ix<nIVals[1];++ix)
        for(iz=0;iz<nIVals[2];++iz)
        {
            ranges[0][ii]=TreePos(tVals[0][ix*2],tVals[1][iy*2],tVals[2][iz*2]);
            ranges[1][ii]=TreePos(tVals[0][ix*2+1],tVals[1][iy*2+1],tVals[2][iz*2+1]);
            ++ii;
        }
        return ranges[0..ii];
    }
    
    enum LoopingFlags:uint{
        Parallel=LoopType.Parallel,
        HalfLoop=Parallel<<1,
    }
    static assert(LoopType.Parallel==1 && LoopType.Sequential == 0,"unexpected values for LoopType");
    
    /// loops on the particles of the given kind within sqrt(cutoff2) of posP1 and in the given ranges
    int neighsRK(int loopType)(SafeV posP1,SafeT cutoff2,TreePos[2][] ranges,KindIdx kind,
        int delegate(ref HPoint!(T)) loopBody,size_t blockSize=1024)
    {
        size_t kMax=kindStarts[ik+1];
        if (kMax>uVal) kMax
        auto kindPos=sortedTreePos[kindStarts[ik]..kindStarts[ik+1]];
        auto kindVals=sortedHPoints[kindStarts[ik]..kindStarts[ik+1]];
        foreach(range;ranges){
            auto lb=lBound(kindPos,range[0]);
            auto ub=uBound(kindPos,range[1]);
            auto rAtt=loopIRange!(loopType)(lb,ub,blockSize1);
            foreach(i;rAtt) {
                auto hAtt=kindVals.ptrI(i);
                SafeV safeP2=toSafeType(hAtt.pos);
                safeP2-=safeP1;
                SafeT d2=safeP2.norm22();
                if (d2<cutoff2){
                    int res=loopBody(hAtt);
                    if (!res) return res;
                }
            }
        }
        return 0;
    }
    /// loops on half the particles of the given kind within sqrt(cutoff2)
    /// avoids self and double counting using B&W counting scheme which balances things well
    /// but is a bit costly..., for efficeincy B&W is done only on the particle
    int neighsRKHalf(int loopType)(PIndex p1,index_type iPos1,SafeV posP1,SafeT cutoff2,TreePos[2][] ranges,KindIdx kind,
        int delegate(ref HPoint!(T)) loopBody,size_t blockSize=1024)
    {
        size_t kMax=kindStarts[ik+1];
        if (kMax>uVal) kMax
        auto kindPos=sortedTreePos[kindStarts[ik]..kindStarts[ik+1]];
        auto kindVals=sortedHPoints[kindStarts[ik]..kindStarts[ik+1]];
        auto lBit=(p1.data&1)!=0;
        foreach(range;ranges){
            auto lb=lBound(kindPos,range[0]);
            auto ub=uBound(kindPos,range[1]);
            auto rAtt=loopIRange!(loopType)(lb,ub,blockSize);
            } else {
                auto rAtt=sLoopIRange(lb,ub);
            }
            foreach(i;rAtt) {
                auto hAtt=kindVals.ptrI(i);
                auto pAtt=hAtt.particle;
                if (pAtt<p1){ // black & white numbering on the particles, switch also on kind parity? would remove "first particle bias"
                    if ((pAtt.data&1)^lBit)!=0) continue;
                } else if (pAtt==p1){
                    if (hAtt.iPos<iPos1){
                        if ((hAtt.iPos+iPos1)&1!=0) continue;
                    } else {
                        if ((hAtt.iPos+iPos1)&1==0) continue;
                    }
                } else {
                    if ((pAtt.data&1)^lBit)==0) continue;
                }
                SafeV safeP2=toSafeType(hAtt.pos);
                safeP2-=safeP1;
                SafeT d2=safeP2.norm22();
                if (d2<cutoff2){
                    int res;
                    res=loopBody(hAtt);
                    if (!res) return res;
                }
            }
        }
        return 0;
    }
    
    /// a structure that loops on the neighbors of the given particle within the given radius
    static struct Looper(int loopingFlags){
        enum { lFlags = loopingFlags }
        enum { lType = ((loopingFlags & LoopingFlags.Parallel)?LoopType.Parallel:LoopType.Sequential) }
        HierarchicalSort hSort;
        size_t blockSize;
        TreePos[2][8] rangesBuf;
        int nRanges;
        HPoint!(T) h1;
        SafeV safeP1;
        /// returns the current ranges
        TreePos[2][] ranges(){
            return rangesBuf[0..nRanges];
        }
        /// sets the current ranges
        void ranges(TreePos[2][]r){
            assert(r.length<=8);
            if (r.ptr !is rangesBuf.ptr){
                rangesBuf[0..r.length]=r;
            }
            nRanges=cast(int)r.length;
        }
        /// constructor
        static Looper opCall(HierarchicalSort hSort,ref HPoint!(T) h1,RVal r,size_t blockSize){
            Looper res;
            res.hSort=hSort;
            res.r=r;
            if (r.valid()) {
                res.ranges=hSort.rangesR(h1.treePos,r,rangesBuf);
            }
            res.h1=h1;
            res.safeP1=toSafeType(h1.pos);
            res.blockSize=blockSize;
            return res;
        }
        /// ditto
        static Looper opCall(HierarchicalSort hSort,PIndex particle1,index_type iPos1,RVal r,size_t blockSize){
            auto h1=sortedHPoints.ptrI(indexOfPos[particle1,iPos1]);
            return opCall(hSort,*h1,r);
        }
        /// loops on the neighbors of the given kind
        int loopKind(KindIdx k,NeighList!(T).LoopBody loopBody){
            if (!r.valid() || ranges == 0) {
                return ellipsoidNeighs!(loopingFlags)(h1.particle,h1.iPos,h1.pos,r,kind,loopBody);
            }
            static if ((loopingFlags & LoopingFlags.HalfLoop)!=0){
                return neighsRKHalf!(lType)(h1.particle,h1.iPos,safeP1,r.cutoff2,ranges,k,
                    delegate int(ref HPoint!(T) h2){
                        NeighPair!(T) pair;
                        pair.p1=h1.particle; // could move out of the loop if sequential
                        pair.i1=h1.iPos; // could move out of the loop if sequential
                        pair.p2=h2.particle;
                        pair.i2=h2.iPos;
                        pair.p1p2=h2.pos-h1.pos;
                        return loopBody(pair);
                    });
            } else {
                return neighsRK!(lType)(safeP1,r.cutoff2,ranges,k,
                    delegate int(ref HPoint!(T) h2){
                        NeighPair!(T) pair;
                        pair.p1=h2.particle;
                        pair.i1=h2.iPos;
                        pair.p2=h1.particle; // could move out of the loop if sequential
                        pair.i2=h1.iPos; // could move out of the loop if sequential
                        pair.p1p2=h1.pos-h2.pos;
                        return loopBody(pair);
                    });
            }
        }
        
        /// loops on all kind of neighbors
        int loopAllKinds(loopBody){
            int res=0;
            foreach (k;kRange.loop!(lType)()){
                if (res!=0) break;
                int rAtt=loopKind(k,loopBody);
                if (rAtt!=0) res=rAtt;
            }
            return res;
        }
    }
    
    /// loops on the neighbors of particle 1 iPos1 with the requested kind that are within r
    int loopOnNeighsPK!(int loopingFlags)(PIndex particle1,index_type iPos1,RVal r,KindIdx kind, NeighList!(T).LoopBody loopBody,
        size_t optSize=defaultOptimalBlockSize)
    {
        auto looper=Looper!(loopingFlags)(this, particle1, iPos1, r, optSize);
        return looper.loopKind(kind,loopBody);
    }
    /// sequential loop on the neighbors of particle 1 iPos1 with the requested kind that are within r
    int sloopOnNeighsPK(PIndex particle1,index_type iPos1,RVal r,KindIdx kind, NeighList!(T).LoopBody loopBody,
        bool avoidDCount=false)
    {
        if (avoidDCount) {
            return loopOnNeighsPK!(LoopingFlags.Sequential|LoopingFlags.HalfLoop)(particle1,iPos1,r,kind,loopBody);
        } else {
            return loopOnNeighsPK!(LoopingFlags.Sequential)(particle1,iPos1,r,kind,loopBody);
        }
    }
    /// parallel loop on the neighbors of particle 1 iPos1 with the requested kind that are within r
    int ploopOnNeighsPK(PIndex particle1,index_type iPos1,RVal r,KindIdx kind, NeighList!(T).LoopBody loopBody,
        bool avoidDCount=false, size_t optSize=defaultOptimalBlockSize)
    {
        if (avoidDCount) {
            return loopOnNeighsPK!(LoopingFlags.Parallel|LoopingFlags.HalfLoop)(particle1,iPos1,r,kind,loopBody);
        } else {
            return loopOnNeighsPK!(LoopingFlags.Parallel)(particle1,iPos1,r,kind,loopBody);
        }
    }

    /// loops on the particle pairs of kinds k1 and k2 that are within r from each other
    int loopKK(int loopingFlags)(KindIdx k1,KindIdx k2,T r,NeighList!(T).LoopBody loopBody,
        size_t optSize=defaultOptimalBlockSize)
    {
        enum { loopType=((loopingFlags&LoopingFlags.Parallel)!=0)?LoopType.Parallel:LoopType.Sequential }
        int resKK=0;
        foreach (ref h1; indexOfPos[k1].loop!(loopType)(optSize)){
            if (resKK != 0) break; // avoid attempt to break out (which sometime has issues), and simply skip iteration?
            auto looper=Looper!(loopingFlags)(this, h1, r, optSize);
            int myRes=looper.loopKind(kind,loopBody);
            if (myRes) resKK=myRes;
        }
        return resKK;
    }
    /// loops on the particle pairs of kinds k1 and k2 that are within r from each other
    int sLoopKK(KindIdx k1,KindIdx k2,T r,NeighList!(T).LoopBody loopBody,bool avoidDCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        if (avoidDCount){
            return loopKK!(LoopingFlags.Sequential|LoopingFlags.HalfLoop)(k1,k2,r,loopBody,optSize);
        } else {
            return loopKK!(LoopingFlags.Sequential)(k1,k2,r,loopBody,optSize);
        }
    }
    /// parallel loop on the particle pairs of kinds k1 and k2 that are within r from each other
    int pLoopKK(KindIdx k1,KindIdx k2,T r,NeighList!(T).LoopBody loopBody,bool avoidDCount=false,
        size_t optSize=defaultOptimalBlockSize)
    {
        if (avoidDCount){
            return loopKK!(LoopingFlags.Parallel|LoopingFlags.HalfLoop)(k1,k2,r,loopBody,optSize);
        } else {
            return loopKK!(LoopingFlags.Parallel)(k1,k2,r,loopBody,optSize);
        }
    }
    
    mixin RefCountMixin!();
}

class HierarchicalSortGen: NeighListGen!(T){
    Real minGridsize;
    
    this(Real minGridsize=0.4){
        this.minGridsize=minGridsize;
    }
    
    /// creates the structure to build neighboring lists on the given configuration 
    NeighList!(Real) createOnConfigReal(Cell!(Real) cell,SegmentedArray!(Vector!(Real,3)) sArr) {
        return new HierarchicalSort!(Real)(cell,sArr,minGridsize);
    }
    /// ditto
    NeighList!(LowP) createOnConfigLowP(Cell!(LowP) cell,SegmentedArray!(Vector!(LowP,3)) sArr) {
        return new HierarchicalSort!(LowP)(cell,sArr,cast(LowP)minGridsize);
    }
    /// ditto
    NeighList!(Real) createOnConfig(Cell!(Real) cell,SegmentedArray!(Vector!(Real,3)) sArr) {
        return new HierarchicalSort!(Real)(cell,sArr,minGridsize);
    }
    /// ditto
    NeighList!(LowP) createOnConfig(Cell!(LowP) cell,SegmentedArray!(Vector!(LowP,3)) sArr) {
        return new HierarchicalSort!(LowP)(cell,sArr,cast(LowP)minGridsize);
    }
}
