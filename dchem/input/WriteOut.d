module dchem.input.WriteOut;
import dchem.Physcon;
import blip.io.BasicIO;
import dchem.sys.ParticleSys;
import dchem.sys.PIndexes;
import dchem.Common;
import dchem.sys.SegmentedArray;

/// writes out an xyz file
void writeXyz(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos,char[] comments){
    // angstrom*...
    auto nAtoms=sysStruct.particles[sysStruct.levels[0]].data.length;
    auto s=dumperP(sink);
    s(nAtoms)("\n");
    foreach(c;comments){
        if (c=='\n' || c=='\r') throw new Exception("comments should have no newline",__FILE__,__LINE__);
    }
    s(comments)("\n");
    foreach(p;sysStruct.externalOrder.lSortedPIndex){
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        s(k.symbol)(" ");
        auto posAtt=pos[p,0];
        s(posAtt.x*angstrom)(" ")(posAtt.y*angstrom)(" ")(posAtt.z*angstrom)("\n");
    }
}

/// writes a turbomole style coordinate file
void writeTurboCoord(T)(CharSink sink,SysStruct sysStruct,SegmentedArray!(Vector!(T,3))pos){
    auto s=dumperP(sink);
    s("$coord\n");
    foreach(p;sysStruct.externalOrder.lSortedPIndex){
        auto k=sysStruct.kinds[cast(size_t)p.kind];
        s(k.name)(" ");
        auto posAtt=pos[p,0];
        s(posAtt.x)(" ")(posAtt.y)(" ")(posAtt.z)(" ")(k.symbol)("\n");
    }
    s("$end\n");
}