/// tests the forces (gradient) using finite differences
module dchem.sampler.TestGrad;
import dchem.Common;
import dchem.sys.ParticleSys;
import dchem.input.RootInput;
import blip.io.BasicIO;
import blip.serialization.Serialization;
import blip.io.Console;
import blip.parallel.smp.PLoopHelpers;
import WriteOut=dchem.input.WriteOut;
import dchem.calculator.FileCalculator;
import dchem.calculator.Calculator;
import dchem.input.WriteOut;
import dchem.sys.DynVars;
import blip.container.GrowableArray;


class TestGrad:Sampler{
    Real dx=0.05;
    Real maxDxFact=2.0;
    Real maxOffsetCFact=0.5;
    Real maxEDiffErr=1.0e-5;
    InputField method;
    mixin myFieldMixin!();
    CharSink _log;
    ubyte[] history;
    ParticleSys!(Real) centralPointReal;
    ParticleSys!(LowP) centralPointLowP;
    DynPVector!(Real,DualDxType) dualForcesReal;// forces (mddx) in the dual space
    DynPVector!(LowP,DualDxType) dualForcesLowP;// forces (mddx) in the dual space

    ParticleSys!(T) centralPointT(T)(){
        static if (is(T==Real)){
            return centralPointReal;
        } else {
            return centralPointLowP;
        }
    }
    
    mixin(serializeSome("dchem.TestGrad",
    `dx: the discretization amount (derivatives are approximated with +dx/2..dx/2)
    maxDxFact: maximum distance between the two evaluation points (needs to be larger than 1)
    maxOffsetCFact: maximum distance from the central point to accept a point evaluation as fraction of dx (needs to be larger than 1)
    maxEDiffErr: maximum error in the energy difference
    method: the method to test`));
    mixin printOut!();
    
    DynPVector!(T,DualDxType) dualForcesT(T)(){
        static if (is(T==Real)){
            return dualForcesReal;
        }
        static if (is(T==LowP)){
            return dualForcesLowP;
        }
    }
    
    /// information about a direction
    class DirInfo{
        TestGrad context;
        Real[2] ens; /// energies evaluated at the different points
        Real[2] dCenter; /// distances fom the central point
        Real coordDx; /// real finite difference distance in the coord space
        Real dualTSDx; /// finite difference distance in the dual tangential space
        Real realTSDx; /// finite difference distance in the real space calculated from the tangential space
        Real cosDir; /// cos of the angle between the requested and actual finite diff direction
        Real gradEDiff; /// expected energy difference using the gradient;
        size_t idir;
        
        this(size_t idir,TestGrad context){
            this.context=context;
            this.idir=idir;
        }
        /// offset from the central point in the coord space
        Real offsetCenter(){
            Real semiP=0.5*(coordDx+dCenter[0]+dCenter[1]);
            Real areaT2=semiP*(semiP-coordDx)*(semiP-dCenter[0])*(semiP-dCenter[1]);
            return 2.0*((areaT2>0)?sqrt(areaT2)/coordDx:0.0);
        }
        /// real energy difference in this direction
        Real eDiff(){
            return ens[0]-ens[1];
        }
        /// error in energy
        Real eError(){
            return eDiff()-gradEDiff;
        }
        /// if the error is too large
        bool eTooLarge(){
            return (coordDx<context.maxDxFact*context.dx &&
                offsetCenter()<context.maxOffsetCFact &&
                eDiff()>context.maxEDiffErr);
        }
        /// writes the header for the data of this direction
        void writeHeader(CharSink sink){
            sink("     dir\tdCenter1\tdCenter2\t offsetC\t coordDx\tdualTSDx\trealTSDx\t  cosDir\t   eDiff\tgradEDiff\teError\tcomments\n");
        }
        /// writes the data regarding this direction
        void writeData(CharSink sink){
            auto s=dumper(sink);
            s(idir,",8"[])("\t")(dCenter[0],",8"[])(dCenter[1],",8"[])("\t")(offsetCenter,",8")("\t");
            s(coordDx,",8"[])("\t")(dualTSDx,",8"[])("\t")(realTSDx,",8"[])("\t")(cosDir,",8")("\t");
            s(eDiff,",8"[])("\t")(gradEDiff,",8"[])("\t")(eError,",8"[])("\t");
            if (eTooLarge) s("ERROR");
            s("\n");
        }
        /// calculates the energies/positions to verify this direction
        void calcE(T)(){
            DynPVector!(T,XType) diff;
            diff=context.centralPointT!(T)().dynVars.dVarStruct.emptyX;
            diff[]=0;
            
            int i=0;
            // calc
            foreach(c;pLoopIter(delegate bool(ref CalculationContext c){
                    if (i<2) {
                        c=method.method.getCalculator(true,context.history);
                        ++i;
                        return true;
                    }
                    return false;
                }))
            {
                auto dualDir=pSysT!(T)(c).dynVars.dVarStruct.emptyDualDx();
                dualDir[]=0;
                T sign=1-2*i;
                dualDir[idir]=sign*0.5*dx;
                auto pippo=&pSysT!(T)(c).addFromDualTSpace!(T);
                pSysT!(T)(c).addFromDualTSpace!(T)(dualDir,pSysT!(T)(c).dynVars.x);
                dualDir.giveBack();
                c.updateEF(true,false);
                synchronized(this){
                    diff.axpby(pSysT!(T)(c).dynVars.x,sign,cast(T)1);
                }
                pSysT!(T)(c).dynVars.x.axpby(centralPointT!(T)().dynVars.x,-1,1);
                dCenter[i]=pSysT!(T)(c).dynVars.x.norm2();
                ens[i]=c.potentialEnergy;
            }
            // check
            coordDx=diff.norm2();
            auto dir=centralPointT!(T)().dynVars.dVarStruct.emptyDx();
            dir[]=0;
            centralPointT!(T)().addToTSpace!(T)(diff,dir);
            gradEDiff=centralPointT!(T)().dotInTSpace(dir,dualForcesT!(T)());

            dir.giveBack();
            auto dualDir=centralPointT!(T)().dynVars.dVarStruct.emptyDualDx();
            centralPointT!(T)().toDualTSpace(dir,dualDir);
            realTSDx=centralPointT!(T)().dotInTSpace(dir,dualDir);
            dualTSDx=dualDir.norm2();
            cosDir=dualDir[idir]/dualTSDx;
        }
    }
    DirInfo[] dirs;
    
    void doAllChecks(T)(){
        sout("pippo doAllChecks\n");
        foreach(i,ref dirI;pLoopArray(dirs,1)){
            sinkTogether(sout,delegate void(CharSink s){ dumper(s)("dir")(i)(" start\n"); });
            dirI=new DirInfo(i,this);
            sinkTogether(sout,delegate void(CharSink s){ dumper(s)("dir")(i)(" preCalc\n"); });
            dirI.calcE!(T)();
            sinkTogether(sout,delegate void(CharSink s){ dumper(s)("dir")(i)(" done\n"); });
        }
        sout("done dir calculations\n");
        if (dirs.length>0){
            dirs[0].writeHeader(sout.call);
        }
        bool failed=false;
        sout("write out data\n");
        foreach(i,dirI;dirs){
            dirI.writeData(sout.call);
            failed=failed||dirI.eTooLarge();
        }
        if (failed){
            sout("Test failed\n");
        } else {
            sout("Test passed\n");
        }
        sout("pippo end doAllChecks\n");
    }
    
    /// does the checks
    void run(){
        
        sout("pippo1\n");
        auto m=method.method(); // method.method;
        if (m is null) throw new Exception("invalid method in field "~myFieldName,__FILE__,__LINE__);
        sout("pippo2\n");
        CalculationContext c=m.getCalculator(true,null);
        sout("pippo3\n");
        size_t nDim;
        sout("pippo4\n");
        mixin(withPSys(`nDim=pSys.dynVars.dVarStruct.dualDxGroup.dataLength;`,"c."));
        sout("pippo5\n");
        dirs=new DirInfo[](nDim);
        sout("pippo6\n");

        sout("pippo7\n");
        c.updateEF(true,true);
        sout("pippo8\n");
        centralPointReal=c.pSysReal();
        if (centralPointReal!is null) centralPointReal=centralPointReal.dup();
        sout("pippo9\n");
        centralPointLowP=c.pSysLowP();
        if (centralPointLowP!is null) centralPointLowP=centralPointLowP.dup();
        sout("pippo10\n");
        auto history=c.storeHistory();
        sout("pippo11\n");
        c.giveBack();
        sout("pippo12\n");
        
        if (centralPointReal!is null){
            doAllChecks!(Real)();
        } else if (centralPointLowP!is null){
            doAllChecks!(LowP)();
        } else {
            throw new Exception("did not get a valid pSys",__FILE__,__LINE__);
        }
        sout("End calculation\n");
    }
    
    /// stops the calculation if possible
    void stop(){
    }
    
    bool verify(CharSink log){
        bool res=true;
        auto s=dumper(log);
        if (method is null || method.typeId!=InputField.TypeId.Method){
            s("Error: method has to be valid and contain a method in field ")(myFieldName)("\n");
            res=false;
        }
        if (!(dx>0)){
            s("Error: dx has to be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        if (!(maxDxFact>1)){
            s("Error: maxDxFact has to be larger than 1 in field ")(myFieldName)("\n");
            res=false;
        }
        if (!(maxOffsetCFact>0)){
            s("Error: maxOffsetCFact has to be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        if (!(maxEDiffErr>0)){
            s("Error: maxEDiffErr has to be larger than 0 in field ")(myFieldName)("\n");
            res=false;
        }
        return res;
    }
}
