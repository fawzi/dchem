#!/bin/bash

calcE=1
calcF=1
self=$0

while [ $# -gt 0 ]
  do
  case $1 in
      --help)
        echo "$self [--help|--e-only|--f-only]"
        ;;
      --e-only)
	calcF=
	;;
      --f-only)
        calcE=
	;;
   esac
   shift
done

if [ -n "$calcE" ] ; then
    cp -f ../../turboCalc/energy ../../turboCalc/ridft.out .
fi
if [ -n "$calcF" ] ; then
    cp -f ../../turboCalc/gradient .
fi
