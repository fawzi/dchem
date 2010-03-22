/// structures to read in external file formats
module dchem.input.ReadIn;
import tango.core.Memory:GC;
import Integer = tango.text.convert.Integer;
import dchem.sys.PIndexes;
import blip.serialization.Serialization;
import blip.serialization.StringSerialize;
import blip.serialization.SerializationMixins;
import blip.text.TextParser;
import dchem.Common;
import tango.text.Regex;
import blip.text.Utils:trim, contains;
import tango.math.Math: min,max;
import blip.io.BasicIO;

/// a residuum particle
struct ResiduumP{
    char[] resName;
    idxType resNr;
    equals_t equals(ResiduumP o){
        return resName==o.resName && resNr==o.resNr;
    }
    hash_t toHash(){
        return getHash(resName,cast(hash_t)resNr);
    }
    int opCmp(ResiduumP o){
        auto r=cmp(resName,o.resName);
        if (r==0) return cmp(resNr,o.resNr);
        return r;
    }
    mixin(serializeSome("","resName|resNr"));
}
/// a chain particle
struct ChainP{
    char[] chainName;
    idxType chainNr;
    equals_t equals(ChainP o){
        return chainName==o.chainName && chainNr==o.chainNr;
    }
    hash_t toHash(){
        return getHash(chainName,cast(hash_t)chainNr);
    }
    int opCmp(ChainP o){
        auto r=cmp(chainName,o.chainName);
        if (r==0) return cmp(chainNr,o.chainNr);
        return r;
    }
    mixin(serializeSome("","chainName|chainNr"));
}

/// represents a read in particle
struct Particle{
    PIndex pIndex;
    PIndex resIndex;
    PIndex chainIndex;
    idxType externalIdx;
    idxType resNr;
    idxType chainNr;
    char[16] _name;
    char[16] _resName;
    char[16] _chainName;
    real[3] pos;
    real charge=0.0; // add occupancy & temperature?
    char[] name(){
        size_t i=_name.length-1;
        while(i!=0 && _name[i]==' ') --i;
        return _name[0..i+1];
    }
    void name(char[]val){
        if (val.length>_name.length){
            throw new Exception("name too long",__FILE__,__LINE__);
        }
        _name[0..val.length]=val;
        for (size_t i=val.length;i<_name.length;++i){
            _name[i]=' ';
        }
    }
    char[] resName(){
        size_t i=name.length-1;
        while(i!=0 && name[i]==' ') --i;
        return name[0..i+1];
    }
    ResiduumP resP(){
        ResiduumP res;
        res.resName=resName;
        res.resNr=resNr;
        return res;
    }
    ChainP chainP(){
        ChainP res;
        res.chainName=chainName;
        res.chainNr=chainNr;
        return res;
    }
    void resName(char[]val){
        if (val.length>_resName.length){
            throw new Exception("resName too long",__FILE__,__LINE__);
        }
        _resName[0..val.length]=val;
        for (size_t i=val.length;i<_resName.length;++i){
            _resName[i]=' ';
        }
    }
    char[] chainName(){
        size_t i=name.length-1;
        while(i!=0 && name[i]==' ') --i;
        return name[0..i+1];
    }
    void chainName(char[]val){
        if (val.length>_chainName.length){
            throw new Exception("chainName too long",__FILE__,__LINE__);
        }
        _chainName[0..val.length]=val;
        for (size_t i=val.length;i<_chainName.length;++i){
            _chainName[i]=' ';
        }
    }
    mixin(serializeSome("",`pIndex|resIndex|chainIndex|externalIdx|resNr|name|resName|chainNr|chainName|pos|charge`));
    mixin printOut!();
}

/// represents a kind
struct Kind{
    char[] name; // name of the kind (should be unique)
    char[] potential; // name of the potential used for the particles
    char[] symbol;
    PIndex nextParticle;
    mixin(serializeSome("",`name|nextParticle|potential|symbol`));
    static Kind opCall(char[] name,PIndex nextParticle){
        Kind res;
        res.name=name;
        res.nextParticle=nextParticle;
        return res;
    }
}

/// represents a read in system of particles
class ReadSystem{
    char[] name; // name of the system
    char[] comments; // comments
    Real[3][3] cell; /// cell info
    Real[3] x0=[0.0,0.0,0.0]; /// lower left pos of the cell
    int[3] periodic; /// periodic directions
    char[] spaceGroup; /// space group
    size_t nParticles; /// number of particles
    Particle[] _particles; // the particles (use particles)
    Kind[] pKinds; /// particle kinds
    Kind[] resKinds; /// residuum kinds
    Kind[] chainKinds; /// chain kinds
    // the following arrays don't use static arrays because it was not supported...
    size_t[char[]] pKindsIdx;     // is linear search more efficient? cache last?
    size_t[char[]] resKindsIdx;   // is linear search more efficient? cache last?
    size_t[char[]] chainKindsIdx; // is linear search more efficient? cache last?
    PIndex[ResiduumP] residui;    /// the residuums particles
    PIndex[ChainP] chains;        /// the chain particles

    static ClassMetaInfo metaI;
    static this(){
        metaI=ClassMetaInfo.createForType!(typeof(this))("dchem.sys.ReadSystem");
        metaI.addFieldOfType!(char[])("name","name of the system");
        metaI.addFieldOfType!(char[])("comments","comments");
        metaI.addFieldOfType!(Real[])("cell","cell matrix");
        metaI.addFieldOfType!(int[])("periodic","which directions are periodic");
        metaI.addFieldOfType!(Real[])("x0","lower left corner of the cell");
        metaI.addFieldOfType!(char[])("spaceGroup","spaceGroup");
        metaI.addFieldOfType!(Particle[])("particles","the particles");
        metaI.addFieldOfType!(Kind[])("pKinds","particle kinds");
        metaI.addFieldOfType!(Kind[])("resKinds","residuum kinds");
        metaI.addFieldOfType!(Kind[])("chainKinds","chain kinds");
    }
    
    ClassMetaInfo getSerializationMetaInfo(){
        return metaI;
    }
    void serial(S)(S s){
        s.field(metaI[0],name);
        s.field(metaI[1],comments);
        Real[] c=(&(cell[0][0]))[0..9];
        s.field(metaI[2],c);
        if (c.ptr !is &(cell[0][0])){
            (&(cell[0][0]))[0..9]=c;
        }
        auto p=periodic[];
        s.field(metaI[3],p);
        if (p.ptr !is periodic.ptr){
            periodic[]=p;
        }
        auto a=x0[];
        s.field(metaI[4],a);
        if (a.ptr !is x0.ptr){
            x0[]=a;
        }
        s.field(metaI[5],spaceGroup);
        auto part=particles;
        s.field(metaI[6],part);
        if (part.ptr !is _particles.ptr){
            particles()[]=part;
        }
        s.field(metaI[7],pKinds);
        s.field(metaI[8],resKinds);
        s.field(metaI[9],chainKinds);
    }
    void serialize(Serializer s){
        serial(s);
    }
    void unserialize(Unserializer s){
        serial(s);
    }
    
    this(size_t initialCapacity=16){
        nParticles=0;
        cell[]=[[1.0,0.0,0.0],
                [0.0,1.0,0.0],
                [0.0,0.0,1.0]];
        periodic[]=0;
        _particles.length=initialCapacity;
        spaceGroup=null;
    }
    
    Particle[] particles(){
        return _particles[0..nParticles];
    }
    
    Particle *addParticle(ref Particle p){
        auto idxAtt=p.name in pKindsIdx;
        if (idxAtt !is null){
            p.pIndex=pKinds[*idxAtt].nextParticle;
            pKinds[*idxAtt].nextParticle +=1;
        } else {
            auto kPos=pKinds.length;
            pKinds~=Kind(p.name.dup,PIndex(cast(KindIdx)kPos,1));
            pKindsIdx[p.name]=kPos;
            p.pIndex=PIndex(cast(KindIdx)kPos,0);
        }
        idxAtt=p.resName in resKindsIdx;
        if (idxAtt !is null){
            auto rp=p.resP;
            auto rIdx=rp in residui;
            if (rIdx is null){
                p.resIndex=resKinds[*idxAtt].nextParticle;
                residui[rp]=p.resIndex;
                pKinds[*idxAtt].nextParticle +=1;
            } else {
                p.resIndex=*rIdx;
            }
        } else {
            auto kPos=resKinds.length;
            resKinds~=Kind(p.resName,PIndex(cast(KindIdx)kPos,1));
            resKindsIdx[p.resName]=kPos;
            p.resIndex=PIndex(cast(KindIdx)kPos,0);
            residui[p.resP]=p.resIndex;
        }
        idxAtt=p.chainName in chainKindsIdx;
        if (idxAtt !is null){
            auto cp=p.chainP;
            auto cIdx=cp in chains;
            if (cIdx is null){
                p.chainIndex=chainKinds[*idxAtt].nextParticle;
                chains[cp]=p.chainIndex;
                pKinds[*idxAtt].nextParticle +=1;
            } else {
                p.chainIndex=*cIdx;
            }
        } else {
            auto kPos=chainKinds.length;
            chainKinds~=Kind(p.chainName,PIndex(cast(KindIdx)kPos,1));
            chainKindsIdx[p.chainName]=kPos;
            p.chainIndex=PIndex(cast(KindIdx)kPos,0);
            chains[p.chainP]=p.chainIndex;
        }
        if (nParticles==_particles.length){
            _particles.length=GC.growLength(nParticles+1,Particle.sizeof)/Particle.sizeof;
        }
        _particles[nParticles]=p;
        ++nParticles;
        return &(_particles[nParticles-1]);
    }
    /// reorders the residuals and chains so that they have the same order as given by resNr, chainNr
    void reorderKinds(){
        // to do
    }
/+    /// deletes the last particle
    void dropLastParticle(){
        if (nParticles==0){
            throw new Exception("dropLastParticle called on empty system",__FILE__,__LINE__);
        }
        --nParticles;
        auto k=&(pKinds[cast(idxType)particles[nParticles].pIndex.kind]);
        if (cast(long)k.nextParticle.particle==1 && (cast(idxType)particles[nParticles].pIndex.kind)+1==pKinds.length){
            pKinds=pKinds[0..pKinds.length-1];
        } else {
            k.nextParticle+=-1;
        }
    }+/
    
    /// clears the particles, but does not dealloc the memory, and the kinds stored, cell,...
    void clearParticles(){
        nParticles=0;
        _particles[]=Particle();
        foreach(i,ref k;pKinds){
            k.nextParticle=PIndex(cast(KindIdx)i,0);
        }
        foreach(i,ref k;resKinds){
            k.nextParticle=PIndex(cast(KindIdx)i,0);
        }
        foreach(i,ref k;chainKinds){
            k.nextParticle=PIndex(cast(KindIdx)i,0);
        }
    }
}

char[] splitEndNr(char[] name,out size_t nr){
    size_t i=name.length;
    while (i!=0){
        --i;
        if (name[i]<'0' || name[i]>'9') break;
    }
    char[] res=name[0..i];
    if (i!=0 && (name[i]=='-'|| name[i]=='_')){
        res=name[0..i-1];
    }
    nr=0;
    char[]nrStr=name[i..$];
    nr=cast(size_t)Integer.parse(nrStr);
    return res;
}

/// reads an xyz frame
ReadSystem readXYZFrame(TextParser!(char) tp,ReadSystem sys=null,bool clearSys=true){
    if (sys is null) sys=new ReadSystem();
    size_t nat;
    try{
        tp(nat);
    } catch(Exception e){
        return null;
    }
    tp.skipLines(1);
    if (clearSys) sys.clearParticles();
    sys.comments=tp.nextLine.dup;
    for (size_t iat=0;iat<nat;++iat){
        Particle p;
        char[] atName,resName,chainName;
        size_t resNr;
        tp.readValue(atName,false);
        p.name=atName;
        tp(p.pos[0])(p.pos[1])(p.pos[2]);
        tp.skipWhitespace();
        if (tp.next(&tp.scanString)){
            resName=splitEndNr(tp.get,resNr);
            p.resNr=cast(idxType)resNr;
            p.resName=resName;
            tp.skipWhitespace();
            if (tp.next(&tp.scanString)){
                chainName=tp.get.dup;
                p.chainName=chainName; // should try to split a number out?
            } else {
                p.chainName="";
            }
        } else {
            p.resName="";
            p.chainName="";
        }
        tp.skipLines(1);
        p.externalIdx=cast(idxType)iat;
        sys.addParticle(p);
    }
    return sys;
}

/// checks if nStr is a number (and only a number), changes D -> E for scientific notation
bool isNumber(char[]nStr){
    size_t i=0;
    bool isNr=false;
    while (i!=nStr.length && nStr[i]==' ')++i;
    if (i!=nStr.length && (nStr[i]=='+'||nStr[i]=='-')) {++i;}
    while (i!=nStr.length && nStr[i]>='0' && nStr[i]<'9'){isNr=true;++i;}
    if (i!=nStr.length && nStr[i]=='.') {++i;}
    while (i!=nStr.length && nStr[i]>='0' && nStr[i]<'9'){isNr=true;++i;}
    if (i!=nStr.length && (nStr[i]=='e' || nStr[i]=='E' || nStr[i]=='d' || nStr[i]=='D')) {
        auto iExp=i;
        ++i;
        isNr=false;
        if (i!=nStr.length && (nStr[i]=='+' || nStr[i]=='-')) ++i;
        while (i!=nStr.length && nStr[i]>='0' && nStr[i]<'9'){isNr=true;++i;}
        while (i!=nStr.length && nStr[i]==' '){ ++i;}
        if (isNr && i==nStr.length){
            if (nStr[iExp]=='d' || nStr[iExp]=='D') nStr[iExp]='E';
            return true;
        }
    }
    while (i!=nStr.length && nStr[i]==' '){ ++i;}
    return isNr && i==nStr.length;
}

struct Splitter(T){
    T[]line;
    dchar[]toSkip;
    int opApply(int delegate(ref size_t ,ref T[])loopBody){
        /+size_t idx=0,pos=0,len=token.length;
        while(pos<len){
            while(line[pos] in toSkip){ ++pos; if (pos ==len) break extLoop; }
            size_t start=pos;
            ++pos;
            while(pos<len && line[pos] !in toSkip){ ++pos; }
            auto res=loopBody(idx,line[start..pos]);
            if (res) return res;
        }+/
        size_t start=0,iword=0;
        int phase=1;
        foreach(idx,c;line){
            if (contains(toSkip,cast(dchar)c)){
                if (phase==2) {
                    auto token=line[start..idx];
                    auto res=loopBody(iword,token);
                    ++iword;
                    if (res) return res;
                    phase=1;
                }
            } else if (phase==1){
                start=idx;
                phase=2;
            }
        }
        if (phase==2){
            auto token=line[start..$];
            auto res=loopBody(iword,token);
            if (res) return res;
        }
    }
    int opApply(int delegate(ref T[])loopBody){
        size_t start=0;
        int phase=1;
        foreach(idx,c;line){
            if (contains(toSkip,cast(dchar)c)){
                if (phase==2) {
                    auto token=line[start..idx];
                    auto res=loopBody(token);
                    if (res) return res;
                    phase=1;
                }
            } else if (phase==1){
                start=idx;
                phase=2;
            }
        }
        if (phase==2){
            auto token=line[start..$];
            auto res=loopBody(token);
            if (res) return res;
        }
    }
    static Splitter opCall(T[] line, dchar[] toSkip){
        Splitter res;
        res.line=line;
        res.toSkip=toSkip;
        return res;
    }
}
Splitter!(T)splitStr(T,U)(T[] line,U[] toSkip){
    return Splitter!(T)(line,cast(dchar[])toSkip);
}
/// reads the header of a car file
/// http://hincklab.uthscsa.edu/~ahinck/html/soft_packs/msi_docs/insight980/formats980/File_Formats_1998.html
ReadSystem readCarHeader(TextParser!(char) tp,ReadSystem sys=null){
    auto line=tp.nextLine();
    auto firstL="!BIOSYM archive 3"; // relax?
    if (line.length<firstL.length || line[0..firstL.length]!=firstL)
        tp.parseError("first line must be '"~firstL~"'",__FILE__,__LINE__);
    line=tp.nextLine();
    auto line2Re = Regex("PBC *= *(ON|OFF|3D|2D|1D)");
    if (! line2Re.test(line))
        tp.parseError("second line must be PBC=ON|OFF|2D'"~firstL~"'",__FILE__,__LINE__);
    switch(line2Re.match(1)){
        case "ON","3D": sys.periodic[]=[1,1,1];
        break;
        case "OFF": sys.periodic[]=[0,0,0];
        break;
        case "2D": sys.periodic[]=[1,1,0];
        break;
        case "1D": sys.periodic[]=[1,0,0];
        break;
        default:
            tp.parseError("error unknown periodicity '"~line2Re.match(1)~"'",__FILE__,__LINE__);
    }
    line=tp.nextLine();
    if (isNumber(line[64..80])){
        sys.comments~="energy="~trim(line[64..80])~"\n";
        sys.name=trim(line[0..64]);
    } else {
        sys.name=trim(line);
    }
    line=tp.nextLine();
    if (line[0..min(line.length,5)]!="!DATE"){
        tp.parseError("expected '!DATE' in the fifth line",__FILE__,__LINE__);
    } else {
        sys.comments~=trim(line)~"\n";
    }

    line=tp.nextLine();
    if (line.length>64){
        sys.name=trim(line[0..64]).dup;
        sys.comments="E="~line[64..$];
    } else {
        sys.name=trim(line).dup;
    }
    if (sys.periodic[2]==1){
        line=tp.nextLine();
        if (line[0..3]!="PBC"){
            tp.parseError("expected PBC line with cell info",__FILE__,__LINE__);
        }
        bool fixedFormat=false;
        if (sys.periodic[0]==1){
            if (line.length>64 && isNumber(line[3..13]) && isNumber(line[13..23])
                && isNumber(line[23..33]) && isNumber(line[33..43]) && isNumber(line[43..53])
                && isNumber(line[53..63]))
            {
                real[6] params;
                int iPos=3;
                size_t ate;
                for (int iparam=0;iparam<6;iparam++){
                    params[iparam]=Float.parse(line[iPos..iPos+10],&ate);
                    if (trim(line[iPos+ate..iPos+10]).length!=0) break;
                    iPos+=10;
                }
                if (iPos==63){
                    auto symmetry=trim(line[63..$]);
                    sys.spaceGroup=symmetry;
                    fixedFormat=true;
                }
            }
            if (!fixedFormat){
                real[6] params;
                foreach (iparam,token;splitStr(line[3..$]," \t\n\r")){
                    if (iparam<6) {
                        params[iparam]=Float.parse(token);
                    } else if (iparam==6){
                        sys.spaceGroup=token;
                        break;
                    }
                }
            }
        } else if (sys.periodic[0]==0 && sys.periodic[1]==1){
            if (line.length>64 && isNumber(line[3..13]) && isNumber(line[13..23])
                && isNumber(line[23..33]))
            {
                real[6] params;
                params[2]=1;
                params[3]=90;
                params[4]=90;
                int iPos=3;
                size_t ate;
                for (int iparam=0;iparam<3;iparam++){
                    params[[0,1,5][iparam]]=Float.parse(line[iPos..iPos+10],&ate);
                    if (trim(line[iPos+ate..iPos+10]).length!=0) break;
                    iPos+=10;
                }
                if (iPos==63){
                    auto symmetry=trim(line[63..$]);
                    sys.spaceGroup=symmetry;
                    fixedFormat=true;
                }
            }
            if (!fixedFormat){
                real[6] params;
                params[2]=1;
                params[3]=90;
                params[4]=90;
                foreach (iparam,token;splitStr(line[3..$]," \t\n\r")){
                    if (iparam<3) {
                        params[[0,1,5][iparam]]=Float.parse(token);
                    } else if (iparam==3){
                        sys.spaceGroup=token;
                        break;
                    }
                }
            }
        } else if (sys.periodic[0]==0 && sys.periodic[1]==0 && sys.periodic[2]==1){
            real[6] params;
            params[1]=1; params[2]=1;
            params[3]=90; params[4]=90; params[5]=90;
            if (line.length>64 && isNumber(line[3..13]))
            {
                int iPos=3;
                size_t ate;
                for (int iparam=0;iparam<1;iparam++){
                    params[iparam]=Float.parse(line[iPos..iPos+10],&ate);
                    if (trim(line[iPos+ate..iPos+10]).length!=0) break;
                    iPos+=10;
                }
                if (iPos==23){
                    auto symmetry=trim(line[63..$]);
                    sys.spaceGroup=symmetry;
                    fixedFormat=true;
                }
            }
            if (!fixedFormat){
                foreach (iparam,token;splitStr(line[3..$]," \t\n\r")){
                    if (iparam<1) {
                        params[iparam]=Float.parse(token);
                    } else if (iparam==1){
                        sys.spaceGroup=token;
                        break;
                    }
                }
            }
        }
    }
    return sys;
}

/// reads an a car frame
ReadSystem readCarFrame(TextParser!(char) tp,ReadSystem sys=null,bool clearSys=true){
    if (sys is null) sys=new ReadSystem();
    size_t iat=0;
    while (1){
// it should be like this, but assuming no section is empty then a simple sequential read should work and support
// also different layouts, change? residue or potential might be just spaces??
//         1-5     atom name        
//         7-20    x Cartesian coordinate of atom      in angstroms  
//         22-35       y Cartesian coordinate of atom      in angstroms  
//         37-50       z Cartesian coordinate of atom      in angstroms  
//         52-55       type of residue containing atom      
//         57-63       residue sequence name       relative to beginning of current molecule, left justified  
//         64-70       potential type of atom      left justified  
//         72-73       element symbol       
//         75-80       partial charge on atom
        // auto line=file.peek(&file.scanLine,true,false); // could be used to perform some checks and maybe use a fixed format
        Particle p;
        char[] atName,resName,chainName;
        size_t resNr;
        tp.readValue(atName,false);
        if (atName=="end") break;
        if (iat==0 && clearSys) sys.clearParticles();
        p.name=atName;
        tp(p.pos[0])(p.pos[1])(p.pos[2]);
        tp.skipWhitespace();
        char[]potString;
        char[] symbol;
        tp(resName);
        resName=splitEndNr(resName,resNr);
        p.resName=resName;
        tp(potString);
        tp(symbol);
        tp(p.charge);
        tp.skipLines(1);
        p.externalIdx=cast(idxType)iat;
        while(true){
            auto kindIdx=p.resName in sys.resKindsIdx;
            if (kindIdx is null) break;
            auto kind=sys.pKinds[*kindIdx];
            if (kind.symbol==symbol && kind.potential==potString) break;
            if (kind.symbol!=symbol){
                throw new Exception("atom name "~p.name~" has two symbols: '"~kind.symbol~"' and '"~symbol~"'",
                    __FILE__,__LINE__);
            }
            // kind.potential != potString
            size_t ic=p.name.length;
            if (ic<p._name.length) {
                p._name[ic]='_';
            } else {
                throw new Exception("could not generate unique names for atom name "~p.name~" with different potentials "~kind.potential~" and "~potString,
                    __FILE__,__LINE__);
            }
            for (size_t ip=0;ip<potString.length;++ip){
                ++ic;
                if (ic>=p.name.length) break; // hope for the best if potential is truncated..., might fail in the next iteration
                p.name[ic]=potString[ip];
            }
        }
        auto pp=sys.addParticle(p);
        auto kind=&(sys.pKinds[cast(size_t)pp.pIndex.kind]);
        if (kind.nextParticle.particle==cast(ParticleIdx)1){
            kind.symbol=symbol;
            kind.potential=potString;
        }
        ++iat;
    }
    if (iat==0) return null;
    return sys;
}

/// read a file in pdb format
ReadSystem readPdb(TextParser!(char) tp,ReadSystem sys=null,bool clearSys=true){
    assert(0,"pdb reading unimplemented");
}
/// reads a frame from the given input
ReadSystem readFrame(MultiReader data,char[] format,long frame=0,ReadSystem readSys=null){
    switch(format){
    case "xyz":
        auto file=new TextParser!(char)(data.readerChar());
        file.newlineIsSpace=false;
        for (long iframe=1;iframe<frame;++iframe){
            size_t nat;
            file(nat);
            file.skipLines(nat+2);
        }
        readSys=readXYZFrame(file,readSys,true);
        if (frame==-1){
            while (true){
                auto readSys1=readXYZFrame(file,readSys,true);
                if (readSys1 is null) break;
                readSys=readSys1;
            }
        }
        break;
    case "car":
        // fully disallow empty frames???
        auto file=new TextParser!(char)(data.readerChar());
        file.newlineIsSpace=false;
        readSys=readCarHeader(file,readSys);
        auto readSys1=readCarFrame(file,readSys,true);
        if (readSys1 !is null) readSys=readSys1;
        for (long iframe=1;iframe<frame;++iframe){
            readSys1=readCarFrame(file,readSys,true);
            if (readSys1 !is null) readSys=readSys1;
        }
        if (frame==-1){
            while (true){
                readSys1=readXYZFrame(file,readSys,true);
                if (readSys1 is null) break;
                readSys=readSys1;
            }
        }
        break;
    case "pdb":
        auto file=new TextParser!(char)(data.readerChar());
        file.newlineIsSpace=false;
        if (frame!=0) throw new Exception("only single frame pdb are supported",__FILE__,__LINE__);
        readSys=readPdb(file,readSys,true);
        break;
    default:
        assert(0,"unexpected fromat "~format);
    }
    return readSys;
}