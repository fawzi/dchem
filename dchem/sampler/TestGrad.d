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
import blip.parallel.smp.Wait;
import blip.math.Math;
import blip.parallel.mpi.MpiModels;

class TestGrad:Sampler{
    Real dx=0.05;
    Real maxDxFact=2.0;
    Real maxOffsetCFact=0.5;
    Real maxEDiffErr=1.0e-5;
    InputField method;
    CharSink _log;
    ubyte[] history;
    ParticleSys!(Real) centralPointReal;
    ParticleSys!(LowP) centralPointLowP;
    DynPVector!(Real,DualDxType) dualForcesReal; // forces (mddx) in the dual space
    DynPVector!(LowP,DualDxType) dualForcesLowP; // forces (mddx) in the dual space

    ParticleSys!(T) centralPointT(T)(){
        static if (is(T==Real)){
            return centralPointReal;
        } else {
            return centralPointLowP;
        }
    }
    
    mixin myFieldMixin!();
    mixin(serializeSome("dchem.TestGrad",`Sampler that tests the gradient using finite differences at the start configuration.`,
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
    static class DirInfo{
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
        /// offset of the finite difference line from the central point in the coord space
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
            return (abs(coordDx)<context.maxDxFact*context.dx &&
                abs(offsetCenter())<context.maxOffsetCFact &&
                abs(eError())>context.maxEDiffErr);
        }
        /// writes the header for the data of this direction
        void writeHeader(CharSink sink){
            sink("     dir\tdCenter1\tdCenter2\t offsetC\t coordDx\tdualTSDx\trealTSDx\t  cosDir\t   eDiff\tgradEDiff\teError\tcomments\n");
        }
        /// writes the data regarding this direction
        void writeData(CharSink sink){
            auto s=dumper(sink);
            s(idir,",8"[])("\t")(dCenter[0],",8"[])("\t")(dCenter[1],",8"[])("\t")(offsetCenter,",8")("\t");
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
            auto diffLock=new RLock();
            
            int i=0;
            bool iterator(ref CalculationContext ctx){
                if (i<2) {
                    ctx=(cast(Method)context.method.contentObj).getCalculator(true,context.history);
                    ++i;
                    return true;
                }
                return false;
            }
            // calc
            foreach(ref c;pLoopIter(&iterator))
            {
                
                auto dualDir=context.centralPointT!(T)().dynVars.dVarStruct.emptyDualDx();
                dualDir[]=0;
                T sign=1-2*(i-1);
                dualDir[idir]=sign*0.5*context.dx;
                auto startP=context.centralPointT!(T)().dup(PSDupLevel.DynProperties|PSDupLevel.DynPNullDx|PSDupLevel.DynPNullMddx);
                startP.addFromDualTSpace!(T)(dualDir,startP.dynVars.x);
                dualDir.giveBack();
                pSysWriterSetT!(T)(c,pSysWriter(startP));
                c.updateEF(true,false);
                {
                    diffLock.lock();
                    scope(exit) diffLock.unlock();
                    diff.opBypax(startP.dynVars.x,sign,cast(T)1);
                }
                startP.dynVars.x.opBypax(context.centralPointT!(T)().dynVars.x,-1,1);
                dCenter[i-1]=startP.dynVars.x.norm2();
                startP.release();
                ens[i-1]=c.potentialEnergy;
                c.giveBack();
            }
            // check
            coordDx=diff.norm2();
            auto dir=context.centralPointT!(T)().dynVars.dVarStruct.emptyDx();
            dir[]=0;
            context.centralPointT!(T)().addToTSpace!(T)(diff,dir);
            gradEDiff= -context.centralPointT!(T)().dotInTSpace(dir,context.dualForcesT!(T)());
            auto dualDir=context.centralPointT!(T)().dynVars.dVarStruct.emptyDualDx();
            context.centralPointT!(T)().toDualTSpace(dir,dualDir);
            realTSDx=context.centralPointT!(T)().dotInTSpace(dir,dualDir);
            realTSDx=((realTSDx>0)?sqrt(realTSDx):0);
            dir.giveBack();
            dualTSDx=dualDir.norm2();
            cosDir=dualDir[idir]/dualTSDx;
        }
    }
    DirInfo[] dirs;
    
    void doAllChecks(T)(){
        sout("evaluated at:")(dynPVectorWriter(centralPointT!(T)().dynVars.x));
        sout("energy:")(centralPointT!(T)().dynVars.potentialEnergy)("\n");
        sout("dualForces:")(dynPVectorWriter(dualForcesT!(T)()));
        foreach(i,ref dirI;pLoopArray(dirs,1)){
            debug(PrintDirCalc) sinkTogether(sout,delegate void(CharSink s){ dumper(s)("dir")(i)(" start\n"); });
            dirI=new DirInfo(i,this);
            dirI.calcE!(T)();
            debug(PrintDirCalc) sinkTogether(sout,delegate void(CharSink s){ dumper(s)("dir")(i)(" done\n"); });
        }
        debug(PrintDirCalc) sout("done dir calculations\n");
        if (dirs.length>0){
            dirs[0].writeHeader(sout.call);
        }
        bool failed=false;
        foreach(i,dirI;dirs){
            dirI.writeData(sout.call);
            failed=failed||dirI.eTooLarge();
        }
        if (failed){
            sout("Test failed\n");
        } else {
            sout("Test passed\n");
        }
    }
    
    /// does the checks
    void run(LinearComm pWorld,CharSink log){
        auto m=cast(Method)method.contentObj(); // method.method;
        if (m is null) throw new Exception("invalid method in field "~myFieldName,__FILE__,__LINE__);
        m.setup(pWorld,log);
        CalculationContext c=m.getCalculator(true,null);
        size_t nDim;
        c.updateEF(true,true);
        auto precision=c.activePrecision();
        switch (precision){
        case Precision.Real:
            centralPointReal=c.refPSysReal();
            centralPointReal=centralPointReal.dup(PSDupLevel.All);
            centralPointReal[]=c.pSysWriterReal();
            nDim=centralPointReal.dynVars.dVarStruct.dualDxGroup.dataLength;
            dualForcesReal=centralPointReal.dynVars.dVarStruct.emptyDualDx();
            centralPointReal.toDualTSpace(centralPointReal.dynVars.mddx,dualForcesReal);
            break;
        case Precision.LowP:
            centralPointLowP=c.refPSysLowP();
            centralPointLowP=centralPointLowP.dup(PSDupLevel.All);
            centralPointLowP[]=c.pSysWriterLowP();
            nDim=centralPointLowP.dynVars.dVarStruct.dualDxGroup.dataLength;
            dualForcesLowP=centralPointLowP.dynVars.dVarStruct.emptyDualDx();
            centralPointLowP.toDualTSpace(centralPointLowP.dynVars.mddx,dualForcesLowP);
            break;
        default:
            assert(0);
        }
        dirs=new DirInfo[](nDim);
        auto history=c.storeHistory();
        c.giveBack();
        
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
        if (method is null || cast(Method)method.contentObj is null){
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
