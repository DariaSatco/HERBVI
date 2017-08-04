DIVON4_SOURCE  = divon4/addbnd.F divon4/bfgs.F divon4/bndopt.F divon4/bndtst.F divon4/bounds.F divon4/bucmve.F       \
                 divon4/bufopt.F divon4/delbnd.F divon4/delete.F divon4/delslv.F divon4/divon.F divon4/dvcopy.F      \
                 divon4/dvdot.F divon4/dvnbkd.F divon4/exmbuc.F divon4/feasmv.F divon4/feqn.F divon4/fun.F divon4/genpnt.F  \
                 divon4/grdcmp.F divon4/intgrl.F divon4/ldlsol.F divon4/locsch.F divon4/modchl.F divon4/mulchk.F divon4/    \
                 divon4/mxstep.F divon4/newptq.F divon4/nmdchl.F divon4/nocut.F divon4/nodaud.F divon4/orthvc.F      \
                 divon4/partn.F divon4/quad.F divon4/quasi.F divon4/ranums.F divon4/recparz.F divon4/rlen.F divon4/settol.F \
                 divon4/shrnk.F divon4/split.F divon4/treaud.F divon4/tredmp.F divon4/tstext.F divon4/usrint.F divon4/usrtrm.F
JIMMY_SOURCE   = jimmy/dfun.F   jimmy/dvnopt.F jimmy/eikfun.F jimmy/eikon.F  jimmy/epfun.F  jimmy/hwmsct.F jimmy/hwnsct.F \
                 jimmy/jimmin.F jimmy/jmarea.F jimmy/jmefin.F jimmy/jmfact.F jimmy/jmgamz.F jimmy/jmgrid.F                \
                 jimmy/jminit.F jimmy/jmkern.F jimmy/jmknif.F jimmy/jmpn.F   jimmy/jmrndm.F jimmy/jmsn.F   jimmy/jmsnfn.F \
                 jimmy/jmxs1.F  jimmy/qcdres.F
CERN_SOURCE    = cernmisc/jmabend.F   cernmisc/jmdbs4.F  cernmisc/jmdbska.F cernmisc/jmmlst.F cernmisc/jmdbesi0.F \
                 cernmisc/jmdbsir3.F  cernmisc/jmmlpt.F  cernmisc/jmucpy.F  cernmisc/jmlnoc.F

#lib_LIBRARIES = libjimmy.a
lib_LTLIBRARIES = libjimmy.la
AM_FFLAGS = -I$(srcdir)/../include -I$(HERWIG6510PATH)/inc -I$(HERWIG6510PATH)/include  -I$(HERWIG6510PATH)/include/herwig6510
libjimmy_la_SOURCES = $(CERN_SOURCE) $(DIVON4_SOURCE) $(JIMMY_SOURCE)

