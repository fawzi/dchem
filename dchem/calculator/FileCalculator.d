module dchem.calculator.FileCalculator;
import blip.io.BasicIO;
import blip.io.StreamConverters;

class TemplateHandler{
    char[][char[]] subs;
    OutWriter[char[]] longSubs;
    VfsFolder sourceDir;
    VfsFolder targetDir;
    size_t maxKeywordLen=128;
    bool shouldWriteReplacementsDict;
    char[] rDictFilename;
    
    /// writes the current replacements dictionary to the given sink
    void writeReplacementsDict(CharSink outF){
        void writeEscaped(char[] str){
            size_t wrote=0;
            foreach(i,c;str){
                switch (c)
                case '\\','\"':
                    outF(str[wrote..i]);
                    outF("\\");
                    wrote=i;
                    break;
                default:
                    break;
            }
            outF(str[wrote..$]);
        }
        outF("{");
        bool cont=false;
        foreach(k,v;subs){
            if (cont) outF(",");
            outF("\n  \"");
            cont=true;
            writeEscaped(k);
            outF("\":\"");
            writeEscaped(v);
            outF("\"");
        }
        foreach(k,v;longSubs){
            if (cont) outF(",");
            outF("\n  \"");
            cont=true;
            writeEscaped(k);
            outF("\":\"");
            v(&writeEscaped);
            outF("\"");
        }
        outF("\n}\n");
    }
    /// evaluates the templates
    void evalTemplates(int level,bool overwriteUnchangedPaths=false){
        if (shouldWriteReplacementsDict){
            auto rDictF=targetDir.file(rDictFilename).create();
            auto outF=strDumper(rDictF.output);
            writeReplacementsDict(outF);
            rDictF.output.close();
        }
        foreach (f;sourceDir.catalog()){
            if (f.name.length>7 && f.name[$-7..$-1]==".templ" &&
                f.name[$-1]>='0' && f.name[$-1]<='9')
            {
                auto lFile=cast(int)(f.name[$-1]-'0');
                auto newName=f.stringOf()[0..$-7];
                auto newName2=;
                if ((lFile==level || (lFile<level && 
                     (!targetDir.file(f.stringOf()[0..$-8]~(cast(char)('0'+level))).exists)) &&
                    !targetDir.file(newName).exists()))
                {
                    auto fIn=f.open.input;
                    scope(exit){
                        fIn.close();
                    }
                    auto newF=targetDir.file(newName).create();
                    auto fOut=newF.output;
                    scope(exit){
                        fOut.close();
                    }
                    makeSubs(inF,strDumper(fOut));
                }
                // 0: all changed, 1: only pos changed, 2: small pos change, 3: extrapolable pos change
            } else if (overwriteUnchangedPaths || level==0 ||
                (!targetDir.file(f.toString()).exists()))
            {
                auto newF=targetDir.file(newName).create(f.open.input);
                auto output=newF.output;
                foreach (l;lIter){
                    makeSubs(l,output);
                }
            }
        }
    }
    
    /// possibly replaces name with a template expansion into outF
    bool maybeReplace(char[] name,CharSink outF){
        auto sub=name in subs;
        if (sub !is null){
            assert((name in longSubs)is null,"double entry");
            outF(*sub);
            return true;
        }
        auto sub2=name in longSubs;
        if (sub2 !is null){
            (*sub2)(outF);
            return true;
        }
        return false;
    }
    
    /// writes the content of inF to outF, performing substitutions on it
    void makeSubs(bool delegate(OutReader) inF,CharSink outF){
        inF(delegate size_t(char[]data, SliceExtent slice,out bool iterate){
            iterate=true;
            for(size_t i=0;i<data.length;++i){
                if (data[i]=='['){
                    if (i==0){
                        for(size_t j=i;j<data.length;++j){
                            if(data[j]=='['){
                                outF(data[0..j]);
                                return j;
                            }
                            if(data[j]==']'){
                                if (!maybeReplace(data[i..(j+1)],outF)) outF(data[i..(j+1)]);
                                return j+1;
                            }
                        }
                        if(data.length<maxKeywordLen){
                            switch (slice){
                            case SliceExtent.Partial:
                                return Eof;
                            case SliceExtent.Maximal:
                                outF(data);
                                throw new IOException("unexpected SliceExtent.Maximal with less than maxKeywordLen size",__FILE__,__LINE__);
                            case ToEnd:
                                iterate=false;
                                outF(data);
                                return data.length;
                            default: assert(0);
                            }
                        }
                        return data.length;
                    } else {
                        outF(data[0..i]);
                        return i;
                    }
                }
            }
            switch (slice){
            case SliceExtent.Partial:
                if (data.length==0)
                    return Eof;
                break;
            case SliceExtent.Maximal:
                if (data.length==0)
                    throw new IOException("unexpected SliceExtent.Maximal with 0 size",__FILE__,__LINE__);
                break;
            case ToEnd:
                iterate=false;
                return data.length;
            default: assert(0);
            }
            outF(data);
            return data.length;
        });
    }
}