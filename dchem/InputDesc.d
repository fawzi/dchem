/// input description
module dchem.InputDesc;
import blip.serialization.Serialization;
import blip.io.BasicIO;
import blip.core.RuntimeTraits;
import blip.core.Array;
import blip.container.Set;
import dchem.input.RootInput;

/// collets and prints out the input description
class InputDesc{
    Set!(ClassMetaInfo) cInfos;

    bool isBasicType(ClassMetaInfo mInfo){
        switch(mInfo.kind){
        case TypeKind.PrimitiveK,TypeKind.ArrayK,TypeKind.DictK,
            TypeKind.AAK,TypeKind.VoidPtr:
            return true;
        default:
            return false;
        }
    }

    this(){
        cInfos=new Set!(ClassMetaInfo)();
    }
    
    void maybeAdd(ClassMetaInfo m){
        if ((!cInfos.contains(m)) && (!isBasicType(m))){
            cInfos.add(m);
            foreach(f;m){
                if (f.metaInfo!is null){
                    maybeAdd(f.metaInfo);
                }
            }
        }
    }
    
    void loadTypes(){
        auto reg=SerializationRegistry();
        synchronized(reg){
            foreach (k,v;reg.name2metaInfos){
                if (k.length>6 && k[0..6]=="input." || (v.ci !is null && implements(v.ci,InputElement.classinfo))){
                    maybeAdd(v);
                }
            }
        }
    }
    
    int opApply(int delegate(ref ClassMetaInfo)loopBody){
        string[] keys=new string[](cInfos.size);
        size_t i=0;
        foreach (k;cInfos){
            keys[i]=k.className;
            ++i;
        }
        assert(i==keys.length);
        sort(keys);
        auto reg=SerializationRegistry();
        foreach(k;keys){
            auto mInfo=reg.name2metaInfos[k];
            auto res=loopBody(mInfo);
            if (res) return res;
        }
        delete keys;
        return 0;
    }
    mixin(descSome("dchem.InputDesc",`cInfos`));
    
    void htmlDesc(CharSink sink){
        auto s=dumper(sink);
s(`
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>Dchem input</title>
            <meta charset="utf-8">
            <link rel="stylesheet" href="http://fawzi.github.com/blip/css/default.css" type="text/css" media="screen">
        </head>
        <body>
            <div class="header"><a href="http://fawzi.github.com/blip"><img src ="http://fawzi.github.com/blip/img/blip100x100.png"></a>
            </div>
    <h1>Dchem input</h1>
    <p>In general the input uses json like syntax:
    <ul>
    <li><b>23</b> for integers</li>
    <li><b>23.5e-4</b>, <b>11.7D+6</b>, <b>17</b>, or <b>6.4</b> for floating point numbers</li>
    <li><b>[1,2,3]</b> for arrays</li>
    <li><b>{key1:value1,key2:value2}</b> for dictionaries</li>
    <li><b>"bla bla"</b> or just <b>strinWithoutSpaces</b> for strings</li>
    <li><b>{class:AClassNameString, field1:value1, field2:value2}</b> for objects.</li>
    </ul>
    </p>
    <p>Unlike json most commas can be dropped. The input consists of a big dictionary, thus it begins with "{", has several key:value elements and ends with "}".</p>
    <p>Here we list of object or structures that are understood by dchem. What you should write as class is the title, and the field names are written in bold. <em>RootInputElement</em> can be written only at the top level of the input, and, excluding Ref, should be written just there.</p>`);
        foreach (cls;this){
            s(`
    <h3><a name="`)(cls.className)(`">`)(cls.className)(`</a></h3>`);
            bool rootInputEl=cls.ci!is null && implements(cls.ci,InputElement.classinfo);
            if (cls.doc.length>0 || rootInputEl){
                s(`
        <p>`);
                if (rootInputEl){
                    s("<em>RootInputElement</em>. ");
                }
                s(cls.doc)("</p>");
            }
            if(cls.nTotFields()>0){
                s(`
    <ul>`);
                foreach(fld;cls){
                    s(`
        <li><a name="`)(cls.className)("-")(fld.name)(`"><b>`)(fld.name)("</b></a>");
                    if (fld.metaInfo){
                        s(" (")(fld.metaInfo.className)(")");
                    }
                    s(":");
                    s(fld.doc);
                    s("</li>");
                }
                s(`
    </ul>`);
            }
        }
        s(`
    </body>
</html>
        `);
    }
}
