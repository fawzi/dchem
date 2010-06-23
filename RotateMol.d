module RotateMol;
import dchem.util.Rotate;
import blip.narray.NArray;
import blip.io.Console;
import Float=tango.text.convert.Float;
import dchem.input.ReadIn;
import tango.io.stream.DataFile;
import blip.text.TextParser;
import blip.io.StreamConverters;


int main(char[][]args){
    if (args.length!=13 && args.length!=14){
        sout("usage: "~args[0]~" x0_1 x0_2 x0_3 axis_1 axis_2 axis_3 v1_1 v1_2 v1_3 v2_1 v2_2 v2_3 [coord.xyz]\n");
        return 1;
    }
    auto x0=zeros!(real)(3),axis=zeros!(real)(3),v1=zeros!(real)(3),v2=zeros!(real)(3);
    int iarg=1;
    for (int i=0;i<3;++i){
        x0[i]=Float.toFloat(args[iarg]);
        ++iarg;
    }
    for (int i=0;i<3;++i){
        axis[i]=Float.toFloat(args[iarg]);
        ++iarg;
    }
    for (int i=0;i<3;++i){
        v1[i]=Float.toFloat(args[iarg]);
        ++iarg;
    }
    for (int i=0;i<3;++i){
        v2[i]=Float.toFloat(args[iarg]);
        ++iarg;
    }
    auto w1=v1-x0;
    auto w2=v2-x0;
    auto n2=norm2(axis);
    if (n2==0) throw new Exception("axis is 0",__FILE__,__LINE__);
    axis/=n2;
    w1.axpby(axis,-dot(w1,axis));
    w2.axpby(axis,-dot(w2,axis));
    auto nw1=norm2(w1);
    auto nw2=norm2(w2);
    if (nw1==0) throw new Exception("v1 is along axis",__FILE__,__LINE__);
    if (nw2==0) throw new Exception("v2 is along axis",__FILE__,__LINE__);
    if (feqrel2(nw1,nw2)<8) {
        serr("WARNING norm of v1 and v2 in rotation frame is quite different:")(nw1)(" vs ")(nw2)("\n");
    }
    w1/=nw1;
    w2/=nw2;
    auto m=eye!(real)(3);
    rotateVV(w1,w2,m);
    serr("res=")(x0.dataPrinter())("+")(m.dataPrinter())("*(x-")(x0.dataPrinter())(")\n");
    if (args.length==14){
        auto file=new TextParser!(char)(toReaderT!(char)((new DataFileInput(args[13])).input));
        file.newlineIsSpace=false;
        auto sys=readXYZFrame(file);
        sout(sys.particles.length)("\n\n");
        foreach (ref p;sys.particles){
            scope pos=a2NA(p.pos);
            pos-=x0;
            pos=dot(m,pos);
            pos+=x0;
            sout(p.name)(" ")(pos[0])(" ")(pos[1])(" ")(pos[2])("\n");
        }
    }
    return 0;
}