#!/bin/bash
eval `/sw/swdist/bin/modulesinit` >& /dev/null
. ~/cp2k/setup.sh >& /dev/null
~/cp2k/exe/Linux-x86-64-intel/cp2k_shell.sopt --port 51000 --id bla >& cp2k.out &
exit 0


