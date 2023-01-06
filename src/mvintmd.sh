#!/usr/bin/env bash
if [ ! -d "./INTMD" ]; then
  mkdir ./INTMD
  printf "Directory \"INTMD\" is created.\n"
else
  printf "Directory \"INTMD\" already exists.\n"
fi
mv f*fit* ./INTMD/
mv crf*fits ./INTMD/
mv acrf*fits ./INTMD/