hsel *fit $I "IMAGETYP=='Bias Frame'" >zero.list
hsel *fit $I "IMAGETYP=='Flat Field'" >flat.list
hsel *fit $I "IMAGETYP=='Light Frame' & EXPTIME > 100 & OBJECT!='fear'" > objall.list
hsel *fit $I "OBJECT=='fear'" > lampall.list
cat flat.list lampall.list objall.list > flatnall.list
cat objall.list lampall.list > specall.list
