#!/usr/bin/env bash
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRCDIR=$BASEDIR/src
python $SRCDIR/makezero_ccdp.py
python $SRCDIR/ccdotz_ccdp.py
python $SRCDIR/makeflat_ccdp.py
python $SRCDIR/makereflat2m.py
python $SRCDIR/divideflat2m.py
python $SRCDIR/removecr_ccdp.py
python $SRCDIR/doapall.py
# python $SRCDIR/identlamp2m.py
python $SRCDIR/reidentlamp2m.py
python $SRCDIR/wavecal2m.py
python $SRCDIR/telluric_base2m.py