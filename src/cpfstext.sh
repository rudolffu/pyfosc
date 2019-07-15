#! /bin/sh
ls YF*fits > templist
cat templist | parallel -j8 python cpfstext.py
if [ ! -d "./raw" ]; then
  mkdir ./raw
fi
mv YF*fits ./raw/
mv templist ./raw/
