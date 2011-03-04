#!/bin/bash

export TURBODIR=/apps/TURBOMOLE62
. $TURBODIR/Config_turbo_env

calcE=1
calcF=1
self=$0

while [ $# -gt 0 ]
  do
  case $1 in
      --help)
        echo "$self [--help|--e-only|--f-only]"
	exit 0
        ;;
      --e-only)
	calcF=
	;;
      --f-only)
        calcE=
	;;
      *)
	echo "unexpected argument $1"
        echo "$self [--help|--e-only|--f-only]"
	exit 1
	;;
   esac
   shift
done

if [ -n "$calcE" ] ; then
    ridft >ridft.out
fi
if [ -n "$calcF" ] ; then
    rdgrad >rdgrad.out
fi
