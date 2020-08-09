#!/usr/bin/env bash
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRCDIR=$BASEDIR/src
# python $SRCDIR/addheaderinfo.py
python $SRCDIR/makezero2m.py
python $SRCDIR/ccdotz.py
python $SRCDIR/makeflat2m.py
python $SRCDIR/makereflat2m.py
python $SRCDIR/divideflat2m.py
python $SRCDIR/removecr2m.py
python $SRCDIR/doapall.py
python $SRCDIR/identlamp2m.py
# python $SRCDIR/identlamp2m.py
python $SRCDIR/wavecal2m.py
python $SRCDIR/telluric_base2m.py