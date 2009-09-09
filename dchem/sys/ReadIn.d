/// structures to read in external file formats
module dchem.sys.ReadIn;
import tango.core.Memory:growLength;
import Integer = tango.text.convert.Integer;

/// represents a read in particle
struct Particle{
    PIndex pIndex;
    PIndex resIndex;
    PIndex chainIndex;
    idxType externalIdx;
    idxType resNr;
    char[16] name;
    char[16] resName;
    char[16] chainName;
    real[3] pos;
}

/// represents a kind
struct Kind{
    char[] resName;
    PIndex nextParticle;
}

/// represents a read in system of particles
class ReadSystem{
    char[] name; // name of the system
    char[] comments; // comments
    real[3][3] cell; /// cell info
    int[3] periodic; /// periodic directions
    char[] spaceGroup; /// space group
    size_t nParticles; /// number of particles
    Particle[] _particles; // the particles (use particles)
    Kind[] pKinds; /// particle kinds
    Kind[] resKinds; /// residuum kinds
    Kind[] chainKinds; /// chain kinds
    size_t[char[16]] pKindsIdx;     // is linear search more efficient? cache last?
    size_t[char[16]] resKindsIdx;   // is linear search more efficient? cache last?
    size_t[char[16]] chainKindsIdx; // is linear search more efficient? cache last?
    
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
    
    void addParticle(ref Particle p){
        auto idxAtt=p.name in pKindsIdx;
        if (idxAtt !is null){
            p.pIndex=pKinds[*idxAtt].nextParticle;
            pKinds[*idxAtt].nextParticle +=1;
        } else {
            auto kPos=pKinds.length;
            pKinds~=Kind(p.name,PIndex(cast(KindIdx)kPos,1));
            pKindsIdx[p.name]=kPos;
            p.pIndex=PIndex(cast(KindIdx)kPos,0);
        }
        idxAtt=p.resName in resKindsIdx;
        if (idxAtt !is null){
            p.resIndex=resKinds[*idxAtt].nextParticle;
            pKinds[*idxAtt].nextParticle +=1;
        } else {
            auto kPos=resKinds.length;
            resKinds~=Kind(p.resName,PIndex(cast(KindIdx)kPos,1));
            resKindsIdx[p.resName]=kPos;
            p.resIndex=PIndex(cast(KindIdx)kPos,0);
        }
        idxAtt=p.chainName in chainKindsIdx;
        if (idxAtt !is null){
            p.chainIndex=chainKinds[*idxAtt].nextParticle;
            pKinds[*idxAtt].nextParticle +=1;
        } else {
            auto kPos=chainKinds.length;
            chainKinds~=Kind(p.chainName,PIndex(cast(KindIdx)kPos,1));
            chainKindsIdx[p.chainName]=kPos;
            p.chainIndex=PIndex(cast(KindIdx)kPos,0);
        }
        if (nParticles==_particles.length){
            _particles.length=growLength(nParticles+1,Particle.sizeof)/Particle.sizeof;
        }
        _particles[nParticles]=p;
        ++nParticles;
    }
    
    /// clears the particles, but does not dealloc the memory, and the kinds stored, cell,...
    void clearParticles(){
        nParticles=0;
        _particles=Particle();
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
ReadSystem readXYZFrame(TextParser!(char) tp,ReadSystem sys=null){
    if (sys is null) sys=new ReadSystem();
    size_t nat;
    tp(nat);
    tp.skipLines(1);
    sys.comments=tp.nextLine().dup();
    for (size_t iat=0;iat<nat;++iat){
        Particle p;
        char[] atName,resName,chainName;
        long resNr;
        real x,y,z;
        tp.readValue(atName,false);
        if(atName.length>p.name.length)
            throw new Exception("particle name too long '"~atName~"'",__FILE__,__LINE__);
        p.name[0..atName.length]=atName;
        p.name[atName.length..$]=' ';
        atName=buf[0..atName.length];
        tp(p.x)(p.y)(p.z);
        tp.skipWhitespace();
        if (tp.next(&tp.scanString)){
            resName=splitEndNr(tp.slice(),resNr);
            if(resName.length>p.resName.length)
                throw new Exception("residuum name too long '"~resName~"'",__FILE__,__LINE__);
            p.resName[0..resName.length]=resName;
            p.resName[resName.length..$]=' ';
            tp.skipWhitespace();
            if (tp.next(&tp.scanString)){
                chainName=tp.slice().dup;
                if(chainName.length>p.chainName.length)
                    throw new Exception("chain name too long '"~chainName~"'",__FILE__,__LINE__);
                p.chainName[0..chainName.length]=chainName;
                p.chainName[chainName.length..$]=' ';
            } else {
                p.chainName[]=' ';
            }
        } else {
            p.resName[]=' ';
            p.chainName[]=' ';
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
    while (i!=nStr.length && nStr[i]==' ')++i;
    if (i!=nStr.length && (nStr[i]=='+'||nStr[i]=='-')) {++i;}
    while (i!=nStr.length && nStr[i]=>'0' && nStr<'9'){isNr=true;++i;}
    if (i!=nStr.length && nStr[i]=='.') {++i;}
    while (i!=nStr.length && nStr[i]=>'0' && nStr<'9'){isNr=true;++i;}
    if (i!=nStr.length && (nStr[i]=='e' || nStr[i]=='E'||nStr[i]=='d' || nStr[i]=='D') {
        auto iExp=i;
        ++i;
        isNr=false;
        if (i!=nStr.length && (nStr[i]=='+' || nStr[i]=='-')) ++i;
        while (i!=nStr.length && nStr[i]=>'0' && nStr<'9'){isNr=true;++i;}
        while (i!=nStr.length && nStr[i]==' '){ ++i;}
        if (isNr && i==nStr.length){
            if (nStr[iExp]=='d' || nStr[iExp]=='D') nStr[iExp]='E';
            return true;
        }
    }
    while (i!=nStr.length && nStr[i]==' '){ ++i;}
    return isNr && i==nStr.length;
}

/// reads the header of a car file
ReadSystem readCarHeader(TextParser!(char) tp,ReadSystem sys=null){
    auto line=tp.nextLine();
    firstL="!BIOSYM archive 3"; // relax?
    if (line.length<firstL.length || line[0..firstL.length]!=firstL)
        throw new ParseException(tp,"first line must be '"~firstL~"'",__FILE__,__LINE__);
    line=tp.nextLine();
    auto line2Re = Regex("PBC *= *(ON|OFF|3D|2D|1D)");
    if (! line2Re.test(line))
        throw new ParseException(tp,"second line must be PBC=ON|OFF|2D'"~firstL~"'",__FILE__,__LINE__);
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
            throw new ParseException(tp,"error unknown periodicity '"~line2Re.match(1)~"'",__FILE__,__LINE__);
    }
    line=tp.nextLine();
    if (line.length>64){
        sys.name=Util.trim(line[0..64]).dup;
        sys.comments="E="~line[64..$];
    } else {
        sys.name=Util.time(line).dup;
    }
    if (sys.periodic[0]==1){
        line=tp.nextLine();
        if (line[0..3]!="PBC"){
            throw new ParseException(tp,"expected PBC line with cell info",__FILE__,__LINE__);
        }
        bool fixedFormat=false;
        if (line.length>64 && isNumber(line[3..13]) && isNumber(line[13..23])
            && isNumber(line[23..33]) && isNumber(line[33..43]) && isNumber(line[43..53])
            && isNumber(line[53..63]))
        {
            real[6] params;
            iPos=3;
            for (int iparam=0;iparam<6;iparam++){
                params[iparam]=Float.parse(line[iPos..iPos+10],&ate);
                if (Util.trim(line[3+ate]).length!=0) break;
                iPos+=10;
            }
            if (iPos==63){
                auto symmetry=Util.trim(line[63..$]);
                sys.comments~="symmetry="~symmetry~"\n";
                if (symmetry.length && symmetry!="none" && symmetry !="c1"){
                    Stdout("WARNING unsupported symmetry ")(symmetry).newline;
                }
                fixedFormat=true;
            }
         }
         if (!fixedFormat)
    }
}

/// reads an a car frame
ReadSystem readCarFrame(TextParser!(char) tp,ReadSystem sys=null){
    if (sys is null) sys=new ReadSystem();
    size_t nat;
    tp(nat);
    tp.skipLines(1);
    auto comments=tp.nextLine();
    for (size_t iat=0;iat<nat;++iat){
        Particle p;
        char[] atName,resName,chainName;
        long resNr;
        real x,y,z;
        tp.readValue(atName,false);
        if(atName.length>p.name.length)
            throw new Exception("particle name too long '"~atName~"'",__FILE__,__LINE__);
        p.name[0..atName.length]=atName;
        p.name[atName.length..$]=' ';
        atName=buf[0..atName.length];
        tp(p.x)(p.y)(p.z);
        tp.skipWhitespace();
        if (tp.next(&tp.scanString)){
            resName=splitEndNr(tp.slice(),resNr);
            if(resName.length>p.resName.length)
                throw new Exception("residuum name too long '"~resName~"'",__FILE__,__LINE__);
            p.resName[0..resName.length]=resName;
            p.resName[resName.length..$]=' ';
            tp.skipWhitespace();
            if (tp.next(&tp.scanString)){
                chainName=tp.slice().dup;
                if(chainName.length>p.chainName.length)
                    throw new Exception("chain name too long '"~chainName~"'",__FILE__,__LINE__);
                p.chainName[0..chainName.length]=chainName;
                p.chainName[chainName.length..$]=' ';
            } else {
                p.chainName[]=' ';
            }
        } else {
            p.resName[]=' ';
            p.chainName[]=' ';
        }
        tp.skipLines(1);
        p.externalIdx=cast(idxType)iat;
        sys.addParticle(p);
    }
    return sys;
}

/// read pdb
