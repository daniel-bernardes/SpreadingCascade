#!/bin/bash

echo P2P NETWORK FILE REQUEST FORMAT CONVERTER > /dev/stderr
echo                                           > /dev/stderr
if [ ! -d "$1" ]; then
  if [ "$1" = "-?" ] || [ "$1" = "--help" ]; then
    echo Converts lists of spreading events { t P C F } to lists > /dev/stderr
    echo of file requests in the format { t C F P1 ... Pn }      > /dev/stderr
  else
    echo Error: \"$1\" is not a valid directory                  > /dev/stderr
  fi	 
  echo                                                             > /dev/stderr
  echo Usage: `basename $0` \<tmp_dir path to support large files\>> /dev/stderr
  echo                                                             > /dev/stderr
  exit;
fi
TMP_SORT="-T "${1}
sort -s -k1,1n -k3,3n -k4,4n -k2,2n $TMP_SORT | awk '
BEGIN { getline; printf "%d %d %d %d",t=$1,client=$3,file=$4,provider=$2; }
{ if (t != $1 || client != $3 || file != $4)
     printf "\n%d %d %d",t=$1,client=$3,file=$4;
  printf " %d",provider=$2; }
END { printf "\n"; }'

