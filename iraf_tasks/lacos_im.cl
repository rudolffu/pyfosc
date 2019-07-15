procedure lacos_im(input,output,outmask)

# cosmic ray rejection script for single images
# by Pieter van Dokkum, April 2001
#
# for HST+WFC data, the following settings seem reasonable:
#   gain=7.
#   readnoise=5.
#   sigclip=4.5
#   sigfrac=0.3
#   objlim=4.
#   niter=4
#
# same settings for HST+STIS, except:
#   gain=1.  (usually)
#   readn=4.4  (usually)
#
# reference: van Dokkum 2001, PASP, 113, 1420
# info at http://www.astro.yale.edu/dokkum/lacosmic/
#
# small update Jan 2012: added "radius=0" to all imreplace commands
#

char input{"",prompt="input image"}
char output{"",prompt="cosmic ray cleaned output image"}
char outmask{"",prompt="output bad pixel map (.pl)"}
real gain{2.,prompt="gain (electrons/ADU) (0=unknown)"}
real readn{6.,prompt="read noise (electrons) (0=unknown)"}
char statsec{"*,*",prompt="section to use for automatic computation of gain"}
real skyval{0.,prompt="sky level that has been subtracted (ADU)"}
real sigclip{4.5,prompt="detection limit for cosmic rays (sigma)"}
real sigfrac{0.5,prompt="fractional detection limit for neighbouring pixels"}
real objlim{1.,prompt="contrast limit between CR and underlying object"}
int niter{4,prompt="maximum number of iterations"}
bool verbose{yes}

begin

 char blk,lapla,deriv2,med5,sub5,sub5abs,med3,med7,imstatout,inputmask
 char noise,sigmap,firstsel,starreject,gfirstsel,finalsel
 char oldoutput,sigima
 char kernel,gkernel
 real midpt,skylev,sig,sigcliplow,usegain
 char l1,l2
 int i
 bool stop
 int previous,npix

 if (!deftask("imcalc")) error(2,"Load package stsdas")
 if (gain<=0 && (!deftask("iterstat"))) error(2,"Supply gain or define task iterstat")

 convolve.bilinear=no
 convolve.radsym=no

 if (verbose) {
  print("")
  print("_______________________________________________________________________________")
  print("")
  print("                 L.A.Cosmic: Laplacian cosmic ray removal")
  print("")
  print("                  P. G. van Dokkum, 2001, PASP 113, 1420")
  print("")
  print("                    Imaging version 1.1 (January 2004)")
  print("_______________________________________________________________________________")
  print("")
  }

 # make temporary files

 blk = mktemp("lacos")
 lapla = mktemp("lacos")
 deriv2 = mktemp("lacos")
 kernel = mktemp("lacos")
 gkernel=mktemp("lacos")
 med3 = mktemp("lacos")
 med5 = mktemp("lacos")
 med7 = mktemp("lacos")
 sub5 = mktemp("lacos")
 sub5abs = mktemp("lacos")
 imstatout = mktemp("lacos")
 noise = mktemp("lacos")
 sigmap = mktemp("lacos")
 firstsel = mktemp("lacos")
 starreject=mktemp("lacos")
 gfirstsel = mktemp("lacos")
 finalsel = mktemp("lacos")
 inputmask = mktemp("lacos")
 oldoutput = mktemp("lacos")
 sigima = mktemp("lacos")

 # create Laplacian kernel

 print("0 -1 0;",>> kernel)
 print("-1 4 -1;",>> kernel)
 print("0 -1 0;",>> kernel)

 # create growth kernel

 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)
 print("1 1 1;",>> gkernel)

 # initialize loop

 usegain = gain

 i=1
 stop=no
 previous=0

 imcopy(input,oldoutput,verb-)
 if (skyval>0) {
  imarith(oldoutput,"+",skyval,oldoutput)
  }
 imcopy(input,outmask,verb-)
 imreplace(outmask,0,upper=INDEF,lower=INDEF,radius=0)

 # start iterations

 while(!stop) {

  if (verbose) {
   print("")
   if (i<10) print("_______________________________ Iteration "//i//" ___________________________________")
   if (i>9) print("_______________________________ Iteration "//i//"___________________________________")
   print("")
   print("")
   }

  # determine gain if unknown; have to assume that observations are sky
  # limited (ie read noise=0)
  # repeat every iteration: gain determination will improve

  if (gain<=0) {

   # determine approximate average sky level

   if (verbose && i==1) print("Trying to determine gain automatically:")
   if (verbose && i>1) print("Improving gain estimate:")
   iterstat(oldoutput//"["//statsec//"]",nsigr=5.,maxiter=10,print-,verb-,
     low=INDEF,upper=INDEF)
   skylev = iterstat.median
   median(oldoutput//"["//statsec//"]",med7,7,7,zlo=INDEF,zhi=INDEF,verb-)
   imarith(oldoutput//"["//statsec//"]","-",med7,sigima)
   imdel(med7)
   imfunc(sigima,sigima,"abs",verb-)
   iterstat(sigima,nsigr=5.,maxiter=10,print-,verb-,low=INDEF,upper=INDEF)
   imdel(sigima)
   sig = iterstat.median * 1.48
   usegain=skylev/sig**2
   if (verbose) {
    print("  Approximate sky level = "//skylev//" ADU")
    print("  Sigma of sky = "//sig)
    print("  Estimated gain = "//usegain)
    print("")
    }
   if (usegain<=0) error(2,"Gain determination failed - provide estimate of gain manually")
   } 

  # take second-order derivative (Laplacian) of input image
  # kernel is convolved with subsampled image, in order to remove negative
  # pattern around high pixels

  if (verbose) {
   print("Convolving image with Laplacian kernel")
   print("")
   }
  blkrep(oldoutput,blk,2,2)
  convolve(blk,lapla,kernel)
  imreplace(lapla,0,upper=0,lower=INDEF,radius=0)
  blkavg(lapla,deriv2,2,2,option="average")

  if (verbose) {
   print("Creating noise model using:")
   print("  gain = "//usegain//" electrons/ADU")
   print("  readnoise = "//readn//" electrons")
   print("")
   }

  # create model of background flux - 5x5 box should exclude all CRs

  median(oldoutput,med5,5,5,zlo=INDEF,zhi=INDEF,verb-)
  imreplace(med5,0.0001,upper=0,lower=INDEF,radius=0)

  # create noise model

  imcalc(med5,noise,"sqrt(im1*"//usegain//" + "//readn//"**2)/"//usegain,verb-)

  # divide Laplacian by noise model

  imarith(deriv2,"/",noise,sigmap)

  # Laplacian of blkreplicated image counts edges twice:

  imarith(sigmap,"/",2.,sigmap)

  # removal of large structure (bright, extended objects)

  imdel(med5)
  median(sigmap,med5,5,5,zlo=INDEF,zhi=INDEF,verb-)
  imarith(sigmap,"-",med5,sigmap)

  # find all candidate cosmic rays
  # this selection includes sharp features such as stars and HII regions

  if (verbose) {
   print("Selecting candidate cosmic rays")
   print("  sigma limit = "//sigclip)
   print("")
   }
  imcopy(sigmap,firstsel,verb-)
  imreplace(firstsel,0,upper=sigclip,lower=INDEF,radius=0)
  imreplace(firstsel,1,lower=0.1,upper=INDEF,radius=0)

  # compare candidate CRs to median filtered image
  # this step rejects bright, compact sources from the initial CR list

  if (verbose) {
   print("Removing suspected compact bright objects (e.g. stars)")
   print("  selecting cosmic rays > "//objlim//" times object flux")
   print("")
   }

  # subtract background and smooth component of objects

  median(oldoutput,med3,3,3,zlo=INDEF,zhi=INDEF,verb-)
  median(med3,med7,7,7,zlo=INDEF,zhi=INDEF,verb-)
  imarith(med3,"-",med7,med3)
  imarith(med3,"/",noise,med3)
  imreplace(med3,0.01,upper=0.01,lower=INDEF,radius=0)

  # compare CR flux to object flux

  imcalc(firstsel//","//sigmap//","//med3,starreject,"(im1*im2)/im3",verb-)

  # discard if CR flux <= objlim * object flux

  imreplace(starreject,0,upper=objlim,lower=INDEF,radius=0)
  imreplace(starreject,1,lower=0.5,upper=INDEF,radius=0)
  imarith(firstsel,"*",starreject,firstsel)

  # grow CRs by one pixel and check in original sigma map

  convolve(firstsel,gfirstsel,gkernel)
  imreplace(gfirstsel,1,lower=0.5,upper=INDEF,radius=0)
  imarith(gfirstsel,"*",sigmap,gfirstsel)
  imreplace(gfirstsel,0,upper=sigclip,lower=INDEF,radius=0)
  imreplace(gfirstsel,1,lower=0.1,upper=INDEF,radius=0)

  # grow CRs by one pixel and lower detection limit

  sigcliplow = sigfrac * sigclip

  if (verbose) {
   print("Finding neighbouring pixels affected by cosmic rays")
   print("  sigma limit = "//sigcliplow)
   print("")
   }

  convolve(gfirstsel,finalsel,gkernel)
  imreplace(finalsel,1,lower=0.5,upper=INDEF,radius=0)
  imarith(finalsel,"*",sigmap,finalsel)
  imreplace(finalsel,0,upper=sigcliplow,lower=INDEF,radius=0)
  imreplace(finalsel,1,lower=0.1,upper=INDEF,radius=0)

  # determine number of CRs found in this iteration

  imdel(gfirstsel)
  imcalc(finalsel//","//outmask,gfirstsel,"(1-im2)*im1",verb-)
  imstat(gfirstsel,fields="npix",lower=0.5,upper=INDEF,for-) | scan(npix)

  # create cleaned output image; use 3x3 median with CRs excluded

  if (verbose) {
   print("Creating output:")
   print("  bad pixel mask: "//outmask)
   print("  cleaned image: "//output)
   print("")
   }
  imdel(med5)
  imarith(outmask,"+",finalsel,outmask)
  imreplace(outmask,1,lower=1,upper=INDEF,radius=0)
  imcalc(outmask,inputmask,"(1.-10000.*im1)",verb-)
  imarith(oldoutput,"*",inputmask,inputmask)
  median(inputmask,med5,5,5,zloreject=-9999,zhi=INDEF,verb-)
  imarith(outmask,"*",med5,med5)
  if (i>1) imdel(output)
  imcalc(oldoutput//","//outmask//","//med5,output,"(1.-im2)*im1+im3",verb-)

  # cleanup and get ready for next iteration

  if (verbose) {
   print("Cleaning up")
   print("")
   }

  imdel(oldoutput)
  imcopy(output,oldoutput,verb-)

  if (verbose) {
   print("")
   print(npix//" cosmic rays found in iteration "//i)
   print("")
   }

  if (npix==0) stop=yes

  i=i+1
  if (i>niter) stop=yes

  # delete temp files

  imdel(blk//","//lapla//","//deriv2//","//med5)
  imdel(med3//","//med7//","//noise//","//sigmap)
  imdel(firstsel//","//starreject//","//gfirstsel)
  imdel(finalsel//","//inputmask)

  }

 if (skyval>0) {
  imarith(output,"-",skyval,output)
  }

 imdel(oldoutput)
 delete(kernel//","//gkernel)

end



