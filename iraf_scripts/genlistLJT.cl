hsel *fits $I "(OBJECT=='Bias_spec')||(OBJECT=='BIas_spec')" >zero.list
hsel *fits $I "OBJECT=='Halogen_slit2.5_G3'" >flat.list
hsel *fits $I "(YGRNM=='Grism 3')" > flatnall.list
hsel *fits $I "(YGRNM=='Grism 3')&(OBJECT!='Halogen_slit2.5_G3')" > specall.list
hsel *fits $I "(YGRNM=='Grism 3')&(OBJECT!='Halogen_slit2.5_G3')&(OBJECT!='Ne+He_slit2.5_G3')&(OBJECT!='He+Ne_slit2.5_G3')" > objall.list
hsel *fits $I "(OBJECT=='Ne+He_slit2.5_G3')||(OBJECT=='He+Ne_slit2.5_G3')" > lampall.list
