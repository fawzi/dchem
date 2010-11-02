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
                for(int isign=-1;isign<2;isign+=2){
                    dAtt[]=cast(Real)0;
                    dAtt[idim]=cast(Real)isign;
                    auto dAtt2=repeat(dAtt,1,1);
                    rotateEiV!(NArray!(Real,1),NArray!(Real,2))(rotOrigin[idir],dir,dAtt2);
                    sout("\ntesting ")(idir)("_")(idim)((isign==-1)?"-":"+")(dAtt.dataPrinter())("\n");
                    foreach(jdir,dir2;dirs){
                        vals[]=dAtt;
                        vals=rotateVEi!(NArray!(Real,1),NArray!(Real,1))(dir2,rotOrigin[idir],vals);
                        unaryOpStr!("*aPtr0=abs(*aPtr0);")(vals);
                        sort(vals.data);
                        sout("+")(jdir)(":")(vals.dataPrinter())("\n");
                    }
                }
            }
        }
    }
    
    this(NArray!(Real,2) dirs,size_t[] rotOrigin, Real expAlpha=1.e3){
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
        auto dist=points-repeat(_pos,points.shape[0],0);
        auto len2=sumAxis(dist*dist,1);
        dist/=repeat(len2,points.shape[1],1);
        auto lVal=empty!(Real)(points.shape);
        foreach(idir,dir;msphereDef.dirs){
            lVal[]=dist;
            rotateVEi(dir,msphereDef.rotOrigin[idir],lVal.T);
            for(index_type idim=0;idim<points.shape[1];++idim){
                auto lightAtt=light[idir][Range(4*idim,4*idim+4)];
                Real lP=lightAtt[0];
                Real lPS=lightAtt[1];
                Real lM=lightAtt[2];
                Real lMS=lightAtt[3];
                index_type[4] indexes=[-1,-2,-3,-4];
                auto lValIdim=lVal[Range(0,-1),idim];
                for(index_type iPoint=0;iPoint<lVal.shape[0];++iPoint){
                    auto lValNow=lValIdim[idim];
                    if (! lValNow<=lPS){
                        if(! lValNow<=lP){
                            lPS=lP;
                            lP=lValNow;
                            indexes[1]=indexes[0];
                            indexes[0]=iPoint;
                            if (isNaN(lValNow)||lP==Real.infinity){
                                lP=Real.infinity;
                                indexes[3]=indexes[2];
                                indexes[2]=iPoint;
                                lMS=lM;
                                lM=-Real.infinity;
                            }
                        } else {
                            lPS=lValNow;
                            indexes[1]=iPoint;
                        }
                    } else if (lValNow<lMS){
                        if (lValNow<lM){
                            lMS=lM;
                            lM=lValNow;
                            indexes[3]=indexes[2];
                            indexes[2]=iPoint;
                            if (lValNow==-Real.infinity){
                                lP=Real.infinity;
                                indexes[1]=indexes[0];
                                indexes[0]=iPoint;
                            }
                        } else {
                            lMS=lValNow;
                            indexes[3]=iPoint;
                        }
                    }
                }
                if (indexes!=[-1,-2,-3,-4]){
                    foreach_reverse(i,idx;indexes){
                        if (idx<0){
                            lightAtt[3-i]=lightAtt[1-idx];
                            auto idim4=idim*4;
                            values[idir,idim4+3-i]=values[idir,idim4+1-idx];
                            origin[idir,idim4+3-i]=origin[idir,idim4+1-idx];
                        } else {
                            lightAtt[3-i]=lVal[idx,idim];
                            auto idim4=idim*4;
                            values[idir,idim4+3-i]=pointValues[idx];
                            origin[idir,idim4+3-i]=dist[idx];
                        }
                    }
                }
            }
        }
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

void testMicrosphere(){
    auto msDef=new MicrosphereDef(a2NAC([//[1.,0.,0.],
        [1./sqrt(2.),1/sqrt(2.),0.],
        [1./sqrt(2.),0.,1/sqrt(2.)],
        [0.,1./sqrt(2.),1./sqrt(2.)],
        [-1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)],
        [1./sqrt(3.),-1./sqrt(3.),1./sqrt(3.)],
        [1./sqrt(3.),1./sqrt(3.),-1./sqrt(3.)],
        [1./sqrt(3.),1./sqrt(3.),1./sqrt(3.)]]),[cast(size_t)0,0,2,0,1,2,0]);
    msDef.evalDirs();
    auto ms=new Microsphere(msDef);
    ms._pos=zeros!(Real)(3);
    ms.addDataPoints(a2NAC([[0.5,0.5,0.5]]),a2NA([1.0]));
    sout("val1:")(ms.value())("\n");
    ms.addDataPoints(a2NAC([[-0.5,-0.5,-0.5]]),a2NA([-1.0]));
    sout("val1:")(ms.value())("\n");
}

void main(){
    testMicrosphere();
}