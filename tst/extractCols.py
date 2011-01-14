#!/bin/env python
import sys

def extract(idx,f,fOut):
    while 1:
        line=f.readline()
        if not line: break
        l=line.split("\t")
        ll=len(l)
        for i in idx:
            if i<ll:
                fOut.write("\t")
                fOut.write(l[i])
            elif not strict:
                fOut.write("\t ")
            else:
                print "index out of bounds for line",repr(line)
        fOut.write("\n")

if __name__=="__main__":
    if len(sys.argv)<4:
        print "extractCols.py fIn fOut index1 [index2 ...]"
        sys.exit(1)
    fIn=file(sys.argv[1])
    fOut=file(sys.argv[2],"w")
    idx=map(int,sys.argv[3:])
    extract(idx,fIn,fOut)
    print "done"
                