## These variables need to exist
prefix=/afs/cern.ch/user/d/dsatco/cern_project/local
exec_prefix=${prefix}
datarootdir=${prefix}/share

## Try to automatically work out the YODA Python path
YODA_PYTHONPATH="/afs/cern.ch/user/d/dsatco/cern_project/local/lib/python2.7/site-packages"
test -n "$YODA_PYTHONPATH" || { (which yoda-config > /dev/null) && YODA_PYTHONPATH=`yoda-config --pythonpath`; }
test -n "$YODA_PYTHONPATH" || echo "yoda-config could not be found: you may need to manually set paths to libYODA and the yoda Python package" 1>&2

export PATH="$exec_prefix/bin:/afs/cern.ch/user/d/dsatco/cern_project/local/lib/../bin:$PATH"
export LD_LIBRARY_PATH="${exec_prefix}/lib:/afs/cern.ch/user/d/dsatco/cern_project/local/lib:/cvmfs/sft.cern.ch/lcg/releases/LCG_85/HepMC/2.06.09/x86_64-slc6-gcc49-opt/lib:/cvmfs/sft.cern.ch/lcg/releases/LCG_85/fastjet/3.2.0/x86_64-slc6-gcc49-opt/lib:/usr/lib:/cvmfs/sft.cern.ch/lcg/releases/LCG_85/MCGenerators/lhapdf/6.1.6/x86_64-slc6-gcc49-opt/lib:$LD_LIBRARY_PATH"
export PYTHONPATH="/afs/cern.ch/user/d/dsatco/cern_project/local/lib/python2.7/site-packages:$YODA_PYTHONPATH:$PYTHONPATH"

export TEXMFHOME="${datarootdir}/Rivet/texmf:$TEXMFHOME"
export HOMETEXMF="${datarootdir}/Rivet/texmf:$HOMETEXMF"
export TEXMFCNF="${datarootdir}/Rivet/texmf/cnf:$TEXMFCNF"
export TEXINPUTS="${datarootdir}/Rivet/texmf/tex//:$TEXINPUTS"
export LATEXINPUTS="${datarootdir}/Rivet/texmf/tex//:$LATEXINPUTS"

if (complete &> /dev/null); then
    test -e "${datarootdir}/Rivet/rivet-completion" && source "${datarootdir}/Rivet/rivet-completion"
fi

unset YODA_PYTHONPATH
