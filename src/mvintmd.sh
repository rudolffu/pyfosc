#!/usr/bin/env bash
if [ ! -d "./INTMD" ]; then
  mkdir ./INTMD
  printf "Directory \"INTMD\" is created.\n"
else
  printf "Directory \"INTMD\" already exists.\n"
fi
mv mask*fits ./INTMD/
mv f*fit* ./INTMD/
mv crf*fits ./INTMD/
mv acrf*fits ./INTMD/
mv 0*fit* ./INTMD/
