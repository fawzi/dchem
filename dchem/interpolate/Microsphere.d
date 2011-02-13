module dchem.interpolate.Microsphere;
import blip.narray.NArray;
import dchem.Common;
import blip.io.Console;
import dchem.util.Rotate;
import blip.math.IEEE;
import blip.math.Math;

class MicrosphereDef{
    NArray!(Real,2) dirs;
    NArray!(Real,2) weights;
    size_t[] rotOrigin;
    Real expAlpha;
    
    index_type nDim(){
        return dirs.shape[1];
    }
    
    void evalDirs(){
        auto dAtt=zeros!(Real)(dirs.shape[1]);
        auto vals=zeros!(Real)(dirs.shape[1]);
        foreach(idir,dir;dirs){
            sout("dir:")(dir.dataPrinter)("\n");
            for(int idim=0;idim<dirs.shape[1];++idim){
                for(int isign=1;isign<2;isign+=2){ // start at -1 if you want to check other dirs, but does not add anything
                    dAtt[]=cast(Real)0;
                    dAtt[idim]=cast(Real)isign;
                    auto dAtt2=repeat(dAtt,1,1);
                    rotateEiV!(NArray!(Real,1),NArray!(Real,2))(rotOrigin[idir],dir,dAtt2);
                    sout("\ntesting ")(idir)("_")(idim)((isign==-1)?"-":"+")(dAtt.dataPrinter())("\n");
                    foreach(jdir,dir2;dirs){
                        vals[]=dAtt;
                        vals=rotateVEi!(NArray!(Real,1),NArray!(Real,1))(dir2,rotOrigin[jdir],vals);
                        unaryOpStr!("*aPtr0=180.0/PI*acos(abs(*aPtr0));")(vals);
                        sort(vals.data);
                        sout("+")(jdir)(":")(vals.dataPrinter())("\n");
                    }
                }
            }
        }
    }
    NArray!(Real,2) dirsArr(){
        auto dirsArr=zeros!(Real)([dirs.shape[0]*dirs.shape[1]*2,dirs.shape[1]]);
        auto dAtt=zeros!(Real)(dirs.shape[1]);
        auto vals=zeros!(Real)(dirs.shape[1]);
        foreach(idir,dir;dirs){
            for(int idim=0;idim<dirs.shape[1];++idim){
                for(int isign=-1;isign<2;isign+=2){
                    dAtt[]=cast(Real)0;
                    dAtt[idim]=cast(Real)isign;
                    auto dAtt2=repeat(dAtt,1,1);
                    rotateEiV!(NArray!(Real,1),NArray!(Real,2))(rotOrigin[idir],dir,dAtt2);
                    auto iidir=idir*dirs.shape[1]*2+idim*2+(isign+1)/2;
                    dirsArr[iidir]=dAtt;
                }
            }
        }
        return dirsArr;
    }
    /// distances using all possible rotation origins
    NArray!(Real,3) allDirDist(){
        auto dMatrix=empty!(Real)([3*dirs.shape[0],3*dirs.shape[0],2]);
        dMatrix[]=90;
        auto dAtt=zeros!(Real)(dirs.shape[1]);
        auto vals=zeros!(Real)(dirs.shape[1]);
        foreach(idir,dir;dirs){
            //for(int isign=-1;isign<2;isign+=2)
            for (size_t iSource=0;iSource<3;++iSource)
            for(int idim=0;idim<dirs.shape[1];++idim){
                dAtt[]=cast(Real)0;
                dAtt[idim]=1.0;
                auto dAtt2=repeat(dAtt,1,1);
                rotateEiV!(NArray!(Real,1),NArray!(Real,2))(iSource,dir/+*(cast(Real)isign)+/,dAtt2);
                foreach(jdir,dir2;dirs)
                //for(int jsign=-1;jsign<2;jsign+=2)
                for (size_t jSource=0;jSource<3;++jSource){
                    if (jdir!=idir){
                        vals[]=dAtt;
                        vals=rotateVEi!(NArray!(Real,1),NArray!(Real,1))(dir2/+*(cast(Real)isign)+/,jSource,vals);
                        unaryOpStr!("*aPtr0=180.0/PI*acos((abs(*aPtr0)>1)?1.0:abs(*aPtr0));")(vals);
                        sort(vals.data);
                        //auto iidir=6*idir+3*(isign+1)/2+iSource;
                        //auto ijdir=6*jdir+3*(jsign+1)/2+jSource;
                        auto iidir=3*idir+iSource;
                        auto ijdir=3*jdir+jSource;
                        if(dMatrix[iidir,ijdir,0]>vals[0]){
                            dMatrix[iidir,ijdir,0]=vals[0];
                            if (dMatrix[iidir,ijdir,1]>vals[1]){
                                dMatrix[iidir,ijdir,1]=vals[1];
                            }
                        } else if (dMatrix[iidir,ijdir,1]>vals[0]){
                            dMatrix[iidir,ijdir,1]=vals[0];
                        }
                    }
                }
            }
        }
        return dMatrix;
    }
    /// distances using the current rotation origins
    NArray!(Real,3) myDirDist(){
        auto dMatrix=empty!(Real)([dirs.shape[0],dirs.shape[0],2]);
        dMatrix[]=90;
        auto dAtt=zeros!(Real)(dirs.shape[1]);
        auto vals=zeros!(Real)(dirs.shape[1]);
        foreach(idir,dir;dirs){
            for(int idim=0;idim<dirs.shape[1];++idim){
                dAtt[]=cast(Real)0;
                dAtt[idim]=1.0;
                auto dAtt2=repeat(dAtt,1,1);
                rotateEiV!(NArray!(Real,1),NArray!(Real,2))(rotOrigin[idir],dir,dAtt2);
                foreach(jdir,dir2;dirs){
                    if (jdir!=idir){
                        vals[]=dAtt;
                        vals=rotateVEi!(NArray!(Real,1),NArray!(Real,1))(dir2,rotOrigin[jdir],vals);
                        unaryOpStr!("*aPtr0=180.0/PI*acos((abs(*aPtr0)>1)?1.0:abs(*aPtr0));")(vals);
                        sort(vals.data);
                        auto iidir=idir;
                        auto ijdir=jdir;
                        if(dMatrix[iidir,ijdir,0]>vals[0]){
                            dMatrix[iidir,ijdir,0]=vals[0];
                            if (dMatrix[iidir,ijdir,1]>vals[1]){
                                dMatrix[iidir,ijdir,1]=vals[1];
                            }
                        } else if (dMatrix[iidir,ijdir,1]>vals[0]){
                            dMatrix[iidir,ijdir,1]=vals[0];
                        }
                    }
                }
            }
        }
        return dMatrix;
    }
    
    this(NArray!(Real,2) dirs,size_t[] rotOrigin, Real expAlpha=1.e3){
        assert(dirs.shape[0]==rotOrigin.length,"dirs and rotOrigin must have the same length");
        this.dirs=dirs;
        this.expAlpha=expAlpha;
        this.rotOrigin=rotOrigin;
        this.weights=ones!(Real)(dirs.shape);
    }
}

class Microsphere{
    MicrosphereDef msphereDef;
    // each dir corresponds to 4 points +,secondary+,-,secondary-
    NArray!(Real,1) _pos;
    NArray!(Real,2) light;
    NArray!(Real,2) values;
    NArray!(Real,3) origin;
    
    this(MicrosphereDef msphereDef){
        this.msphereDef=msphereDef;
        auto nDim=msphereDef.nDim;
        auto nDir=msphereDef.dirs.shape[0];
        _pos=zeros!(Real)(nDim);
        light=zeros!(Real)([nDir,4*nDim]);
        values=zeros!(Real)([nDir,4*nDim]);
        origin=zeros!(Real)([nDir,4*nDim,nDim]); // avoid storing?
    }
    
    void setPos(NArray!(Real,1)pos){
        this._pos=pos;
        light[]=cast(Real)0;
        // no really needed
        values[]=Real.init;
        origin[]=Real.init;
    }
    
    void addDataPoints(NArray!(Real,2)points,NArray!(Real,1) pointValues){
        sout("pippo1\n");
        auto dist=points-repeat(_pos,points.shape[0],0);
        auto len2=sumAxis(dist*dist,1);
        dist/=repeat(len2,points.shape[1],1);
        auto lVal=empty!(Real)(points.shape);
        sout("pippo2\n");
        foreach(idir,dir;msphereDef.dirs){
            sout("pippo3\n");
            lVal[]=dist;
            rotateVEi(dir,msphereDef.rotOrigin[idir],lVal.T);
            sout("pippo4\n");
            for(index_type idim=0;idim<points.shape[1];++idim){
                sout("pippo5\n");
                auto lightAtt=light[idir][Range(4*idim,4*idim+4)];
                Real lP=lightAtt[0];
                Real lPS=lightAtt[1];
                Real lM=lightAtt[2];
                Real lMS=lightAtt[3];
                sout("pippo6\n");
                index_type[4] indexes=[-1,-2,-3,-4];
                auto lValIdim=lVal[Range(0,-1),idim];
                sout("pippo7\n");
                for(index_type iPoint=0;iPoint<lVal.shape[0];++iPoint){
                    sout("pippo8\n");
                    auto lValNow=lValIdim[idim];
                    if (! lValNow<=lPS){
                        if(! lValNow<=lP){
                            sout("pippo9\n");
                            lPS=lP;
                            lP=lValNow;
                            indexes[1]=indexes[0];
                            indexes[0]=iPoint;
                            sout("pippo10\n");
                            if (isNaN(lValNow)||lP==Real.infinity){
                                sout("pippo11\n");
                                lP=Real.infinity;
                                indexes[3]=indexes[2];
                                indexes[2]=iPoint;
                                lMS=lM;
                                lM=-Real.infinity;
                                sout("pippo12\n");
                            }
                        } else {
                            sout("pippo13\n");
                            lPS=lValNow;
                            indexes[1]=iPoint;
                            sout("pippo14\n");
                        }
                    } else if (lValNow<lMS){
                        if (lValNow<lM){
                            sout("pippo15\n");
                            lMS=lM;
                            lM=lValNow;
                            indexes[3]=indexes[2];
                            indexes[2]=iPoint;
                            sout("pippo16\n");
                            if (lValNow==-Real.infinity){
                                sout("pippo17\n");
                                lP=Real.infinity;
                                indexes[1]=indexes[0];
                                indexes[0]=iPoint;
                                sout("pippo18\n");
                            }
                        } else {
                            sout("pippo19\n");
                            lMS=lValNow;
                            indexes[3]=iPoint;
                            sout("pippo20\n");
                        }
                    }
                }
                sout("pippo21\n");
                if (indexes!=[-1,-2,-3,-4]){
                    foreach_reverse(i,idx;indexes){
                        if (idx<0){
                            sout("pippo22 ")(3-i)(" ")(1-idx)(" lightAtt:")(lightAtt)(" indexes:")(indexes)("\n");
                            lightAtt[3-i]=lightAtt[-1-idx];
                            sout("pippo22a\n");
                            auto idim4=idim*4;
                            sout("pippo22b\n");
                            values[idir,idim4+3-i]=values[idir,idim4+1-idx];
                            sout("pippo22c\n");
                            origin[idir,idim4+3-i]=origin[idir,idim4+1-idx];
                            sout("pippo22d\n");
                        } else {
                            sout("pippo23\n");
                            lightAtt[3-i]=lVal[idx,idim];
                            auto idim4=idim*4;
                            values[idir,idim4+3-i]=pointValues[idx];
                            origin[idir,idim4+3-i]=dist[idx];
                            sout("pippo24\n");
                        }
                    }
                }
            }
        }
        sout("pippo25\n");
    }
    
    Real value(){
        auto nDim=msphereDef.nDim;
        Real wTot=0;
        Real val=0;
        foreach(idir,lVal;light){
            scope vAtt=values[idir];
            for (index_type idim=0;idim<nDim;++idim){
                auto idim0=idim*4;
                auto lAtt=lVal[idim0];
                if (lAtt!=0){
                    wTot+=lAtt;
                    val+=lAtt*values[idir,idim0];
                    auto lSecond=lVal[idim0+1];
                    if (lSecond!=0){
                        lSecond*=exp(msphereDef.expAlpha*(1.0-lAtt/lSecond));
                        if (lSecond!=0){
                            wTot+=lSecond;
                            val+=lSecond*vAtt[idim0+1];
                        }
                    }
                }
                lAtt=lVal[idim0+2];
                if (lAtt!=0){
                    wTot-=lAtt;
                    val-=lAtt*vAtt[idim0+2];
                    auto lSecond=lVal[idim0+3];
                    if (lSecond!=0){
                        lSecond*=exp(msphereDef.expAlpha*(1.0-lAtt/lSecond));
                        if (lSecond!=0){
                            wTot-=lSecond;
                            val-=lSecond*vAtt[idim0+3];
                        }
                    }
                }
            }
        }
        return val/wTot;
    }
}

NArray!(Real,2) default3SphereDirs;
size_t[] default3SphereROrigin=[cast(size_t)0,0,0,0,0];

static this(){
    default3SphereDirs=a2NAC([[1.,0.,0.],
        [-1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)],
        [1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)],
        [1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)],
        [1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)]]);
}

T binom(T)(T n,T k){
    T res=1;
    for (T i=k+1;i<=n;++i) res*=i;
    for (T i=n-k;i>0;--i) res/=i;
    return res;
}

/// this generates the +/-1,..,+/-1 directions to rotate to to sample the hypersphere.
/// there are 2^(dim-1) unique directions, this this approach is doable only for a small
/// number of directions.
/// I suspect that when dim is odd last rotations are already done, but ti should be checked
/// not yet tested
NArray!(T,2) generate1Dirs(T)(size_t dim){
    if (dim==0 || dim==1)
        return zeros!(T)([0,dim]);
    auto nDirs=powI(2,dim-2);
    if (dim%2==0){
        nDirs-=binom(dim,dim/2);
    }
    auto res=empty!(T)([nDirs,dim]);
    size_t nNeg=0;
    auto iStack=new size_t[](dim);
    dirAtt=empty!(T)(dim);
    T val=cast(T)1./sqrt(dim);
    res[]=val;
    size_t iidir=1;
    for(size_t nNeg=1;nNeg<(dim-1)/2;++nNeg){
        for (size_t i=0;i<nNeg;++i){
            iStack[i]=i;
        }
        while(true){
            auto iNow=nNeg-1;
            ++iStack[iNow];
            while(iNow>0 && iStack[iNow]>dim-(nNeg-iNow-1)){
                ++iStack[--iNow];
            }
            if (iNow==0 && iStack[0]==dim-nNeg) break;
            for (size_t i=iNow+1;i<nNeg;++i){
                iStack[i]=iStack[i-1]+1;
            }
            foreach(i;iStack){
                dirAtt[iidir,i]=-val;
            }
            ++iidir;
        }
    }
/+    if (dim%2==0){
    // in 2d all these directions are already covered, I suspect it is the case even in more dims 
    // but I should check it
        auto nNeg=(dim-1)/2;
        for (size_t i=0;i<nNeg;++i){
            iStack[i]=i+1;
        }
        while(true){
            auto iNow=nNeg-1;
            ++iStack[iNow];
            while(iNow>0 && iStack[iNow]>dim-(nNeg-iNow)){
                ++iStack[--iNow];
            }
            if (iNow==0 && iStack[0]==dim-nNeg) break;
            for (size_t i=iNow+1;i<nNeg;++i){
                iStack[i]=iStack[i-1]+1;
            }
            dirAtt[iidir,0]=-val;
            foreach(i;iStack){
                dirAtt[iidir,i]=-val;
            }
            ++iidir;
        }
    }+/
    assert(iidir==res.shape[0]);
    return res;
}

/+allSymm=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],
    [1./sqrt(2.),1/sqrt(2.),0.],
    [1./sqrt(2.),0.,1/sqrt(2.)],
    [0.,1./sqrt(2.),1./sqrt(2.)],
    [-1./sqrt(2.),1/sqrt(2.),0.],
    [-1./sqrt(2.),0.,1/sqrt(2.)],
    [0.,-1./sqrt(2.),1./sqrt(2.)],
    [-1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)],
    [1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)],
    [1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)],
    [1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)]];+/

void testMicrosphere(){
    auto msDef=new MicrosphereDef(default3SphereDirs,default3SphereROrigin);
    auto dists=msDef.myDirDist();
    auto d1=dists[Range(0,-1),Range(0,-1),0];
    auto d2=dists[Range(0,-1),Range(0,-1),1];
    sout("dists\n")(d1.dataPrinter("F0,2",80))("\nsecondDist:\n")(d1.dataPrinter("F0,2",80))("\n");
    //msDef.evalDirs();
    auto dirs=msDef.dirsArr();
    sout("dirs:\n")(dirs.dataPrinter("F6,10"))("\n");
    auto ms=new Microsphere(msDef);
    ms._pos=zeros!(Real)(3);
    ms.addDataPoints(a2NAC([[0.5,0.5,0.5]]),a2NA([1.0]));
    sout("val1:")(ms.value())("\n");
    ms.addDataPoints(a2NAC([[-0.5,-0.5,-0.5]]),a2NA([-1.0]));
    sout("val1:")(ms.value())("\n");
}

version(mkMain){
void main(){
    testMicrosphere();
}
}
