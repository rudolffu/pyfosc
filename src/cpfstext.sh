#! /bin/sh
ls YF*fits > templist
{
  command -v parallel && cat templist | parallel -j8 cpfstext.py
} || {
for file in YF*fits 
do
    cpfstext.py $file
done
}
if [ ! -d "./raw" ]; then
  mkdir ./raw
fi
mv YF*fits ./raw/
mv templist ./raw/
