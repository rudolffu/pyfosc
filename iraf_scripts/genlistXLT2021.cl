hsel *fit $I "OBSTYPE=='BIAS'" >zero.list
hsel *fit $I "OBSTYPE=='SPECLFLAT'" >flat.list
hsel *fit $I "OBSTYPE=='SPECLTARGET'" > objall.list
hsel *fit $I "OBSTYPE=='SPECLWLLIGHT'" > lampall.list
cat flat.list lampall.list objall.list > flatnall.list
cat objall.list lampall.list > specall.list
