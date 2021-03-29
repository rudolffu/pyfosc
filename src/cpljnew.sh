#! /bin/sh
# if [ ! -d "./fz_bak" ]; then
#   mkdir ./fz_bak
# fi
# ls *.fz | parallel funpack {}
# mv *fits.fz ./fz_bak/

ls lj*fits > templist
{
  command -v parallel && cat templist | parallel -j8 cpljnew.py
} || {
for file in lj*fits 
do
    cpljnew.py $file
done
}
if [ ! -d "./raw" ]; then
  mkdir ./raw
fi
mv lj*fits ./raw/
mv templist ./raw/
