## Process this file with automake to produce Makefile.in

SUBDIRS = include/ngspice include/cppduals/duals misc maths frontend spicelib
DIST_SUBDIRS = include/ngspice include/cppduals/duals misc maths frontend spicelib xspice ciderlib osdi

if XSPICE_WANTED
SUBDIRS += xspice
endif

if OSDI_WANTED
SUBDIRS += osdi
endif

if CIDER_WANTED
SUBDIRS += ciderlib
endif

if !TCL_MODULE
if !SHARED_MODULE
bin_PROGRAMS =  ngspice

if OLDAPPS
if !WINGUI
bin_PROGRAMS += ngnutmeg ngsconvert ngproc2mod ngmultidec ngmakeidx
helpdatadir = $(pkgdatadir)/helpdir
helpdata_DATA = ngspice.txt ngspice.idx
if !NO_HELP
bin_PROGRAMS += nghelp
endif !NO_HELP
endif !WINGUI
endif OLDAPPS
endif !SHARED_MODULE
endif !TCL_MODULE

EXTRA_DIST = ngspice.txt setplot spectrum \
	     devload devaxis ciderinit \
	     hist_info.c hist_info.h winmain.c winmain.h \
	     tclspice.c tclspice.map pkgIndex.tcl.in spinit.in tclspinit.in \
	     main.c sharedspice.c ngicon.rc ngspice.ico


initdatadir = $(pkgdatadir)/scripts

initdata_DATA = spinit setplot spectrum
if CIDER_WANTED
initdata_DATA += devload devaxis ciderinit
endif

DYNAMIC_DEVICELIBS =					\
	spicelib/devices/asrc/libasrc.la		\
	spicelib/devices/bjt/libbjt.la			\
	spicelib/devices/bsim1/libbsim1.la		\
	spicelib/devices/bsim2/libbsim2.la		\
	spicelib/devices/bsim3/libbsim3.la		\
	spicelib/devices/bsim3v0/libbsim3v0.la		\
	spicelib/devices/bsim3v1/libbsim3v1.la		\
	spicelib/devices/bsim3v32/libbsim3v32.la	\
	spicelib/devices/bsim4/libbsim4.la		\
	spicelib/devices/bsim4v5/libbsim4v5.la		\
	spicelib/devices/bsim4v6/libbsim4v6.la		\
	spicelib/devices/bsim4v7/libbsim4v7.la		\
	spicelib/devices/cap/libcap.la			\
	spicelib/devices/bsim3soi_pd/libbsim3soipd.la	\
	spicelib/devices/bsim3soi_fd/libbsim3soifd.la	\
	spicelib/devices/bsim3soi_dd/libbsim3soidd.la	\
	spicelib/devices/bsimsoi/libbsim4soi.la		\
	spicelib/devices/cccs/libcccs.la		\
	spicelib/devices/ccvs/libccvs.la		\
	spicelib/devices/cpl/libcpl.la			\
	spicelib/devices/csw/libcsw.la			\
	spicelib/devices/dio/libdio.la			\
	spicelib/devices/ind/libind.la			\
	spicelib/devices/isrc/libisrc.la		\
	spicelib/devices/hfet1/libhfet.la		\
	spicelib/devices/hfet2/libhfet2.la		\
	spicelib/devices/hicum2/libhicum2.la		\
	spicelib/devices/hisim2/libhisim2.la		\
	spicelib/devices/hisimhv1/libhisimhv1.la	\
	spicelib/devices/hisimhv2/libhisimhv2.la	\
	spicelib/devices/jfet/libjfet.la		\
	spicelib/devices/jfet2/libjfet2.la		\
	spicelib/devices/ltra/libltra.la		\
	spicelib/devices/mes/libmes.la			\
	spicelib/devices/mesa/libmesa.la		\
	spicelib/devices/mos1/libmos1.la		\
	spicelib/devices/mos2/libmos2.la		\
	spicelib/devices/mos3/libmos3.la		\
	spicelib/devices/mos6/libmos6.la		\
	spicelib/devices/mos9/libmos9.la		\
	spicelib/devices/res/libres.la			\
	spicelib/devices/soi3/libsoi3.la		\
	spicelib/devices/sw/libsw.la			\
	spicelib/devices/txl/libtxl.la			\
	spicelib/devices/tra/libtra.la			\
	spicelib/devices/urc/liburc.la			\
	spicelib/devices/vbic/libvbic.la		\
	spicelib/devices/vccs/libvccs.la		\
	spicelib/devices/vcvs/libvcvs.la		\
	spicelib/devices/vdmos/libvdmos.la		\
	spicelib/devices/vsrc/libvsrc.la




## Build ngspice first:

## compile the icon:
if WINRESOURCE
ngicon.o: ngicon.rc
	windres $(srcdir)/ngicon.rc -O coff -o ngicon.o
endif

ngspice_SOURCES = \
	main.c \
	conf.c \
	conf.h \
	ngspice.c

ngspice_CPPFLAGS = $(AM_CPPFLAGS) -DSIMULATOR

if WINGUI
ngspice_LDFLAGS = -municode $(AM_LDFLAGS)
ngspice_SOURCES += winmain.c hist_info.c
else
ngspice_LDFLAGS = 
endif

ngspice_LDADD = \
	frontend/libfte.la

if WINGUI
ngspice_LDADD += \
	frontend/wdisp/libwindisp.la
endif

ngspice_LDADD += \
	frontend/plotting/libplotting.la

ngspice_LDADD += \
	spicelib/devices/dev.lo \
	$(DYNAMIC_DEVICELIBS)

if NDEV_WANTED
ngspice_LDADD += spicelib/devices/ndev/libndev.la
endif

if NUMDEV_WANTED
ngspice_LDADD += \
	spicelib/devices/nbjt/libnbjt.la   \
	spicelib/devices/nbjt2/libnbjt2.la \
	spicelib/devices/numd/libnumd.la   \
	spicelib/devices/numd2/libnumd2.la \
	spicelib/devices/numos/libnumos.la
endif

ngspice_LDADD += \
	spicelib/analysis/libckt.la \
	spicelib/devices/libdev.la

if XSPICE_WANTED
ngspice_LDADD += \
	xspice/cm/libcmxsp.la \
	xspice/mif/libmifxsp.la \
	xspice/evt/libevtxsp.la \
	xspice/enh/libenhxsp.la \
	xspice/ipc/libipcxsp.la \
	xspice/idn/libidnxsp.la 
endif
ngspice_LDADD += $(XSPICEDLLIBS)

ngspice_LDADD += \
	frontend/parser/libparser.la \
	frontend/numparam/libnumparam.la \
	frontend/trannoise/libtrannoise.la \
	spicelib/parser/libinp.la

if CIDER_WANTED
ngspice_LDADD += \
	ciderlib/twod/libcidertwod.la     \
	ciderlib/oned/libcideroned.la     \
	ciderlib/input/libciderinput.la   \
	ciderlib/support/libcidersuprt.la
endif

if OSDI_WANTED
ngspice_LDADD += osdi/libosdi.la 
endif

ngspice_LDADD += \
	maths/deriv/libderiv.la \
	maths/cmaths/libcmaths.la \
	maths/misc/libmathmisc.la \
	maths/fft/libmathfft.la \
	maths/poly/libpoly.la \
	maths/ni/libni.la

if KLU_WANTED
ngspice_LDADD += maths/KLU/libKLU.la
endif

ngspice_LDADD += \
	maths/sparse/libsparse.la \
	misc/libmisc.la

if SP_WANTED
ngspice_LDADD += \
    maths/dense/libdense.la
endif

if !WINGUI
if !NO_X
ngspice_LDADD += frontend/help/libhlp.la
endif
endif

if WINRESOURCE
ngspice_SOURCES += ngicon.rc
endif

.rc.lo:
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --tag=RC --mode=compile $(RC) $(RCFLAGS) $< -o $@

if OLDAPPS
## nutmeg:

ngnutmeg_SOURCES = \
	main.c \
	conf.c \
	conf.h \
	ngnutmeg.c

ngnutmeg_CPPFLAGS = $(AM_CPPFLAGS)

if WINGUI
ngnutmeg_SOURCES += winmain.c hist_info.c
endif

ngnutmeg_LDADD = \
	frontend/libfte.la

if WINGUI
ngnutmeg_LDADD += \
	frontend/wdisp/libwindisp.la
endif

ngnutmeg_LDADD += \
	frontend/plotting/libplotting.la \
	frontend/parser/libparser.la \
	frontend/numparam/libnumparam.la \
	frontend/trannoise/libtrannoise.la \
	maths/cmaths/libcmaths.la \
	maths/misc/libmathmisc.la \
	maths/fft/libmathfft.la \
	maths/poly/libpoly.la \
	misc/libmisc.la \
	spicelib/parser/libinp.la


## These programs are not required on Windows
if !WINGUI
if !SHWIN

if !NO_X
ngnutmeg_LDADD += \
	frontend/help/libhlp.la
endif !NO_X

## help:

nghelp_SOURCES = conf.c nghelp.c


if !NO_HELP
# We keep the libraries to a minimum here and link against terminal.o
# rather than the full front-end library libfte.a to avoid link errors that
# that would otherwise occur (thanks to Andreas Unger for this fix).
nghelp_LDADD = \
	frontend/terminal.lo \
	misc/libmisc.la
if !NO_X
nghelp_LDADD += \
	frontend/help/libhlp.la
endif !NO_X
endif !NO_HELP


## sconvert:

ngsconvert_SOURCES =  ngsconvert.c
ngsconvert_LDADD = \
	frontend/dvec.lo \
    frontend/rawfile.lo \
    frontend/typesdef.lo \
    frontend/variable.lo \
    frontend/dimens.lo \
    frontend/quote.lo \
    frontend/parser/lexical.lo \
    frontend/parser/numparse.lo \
    misc/string.lo \
    misc/wlist.lo \
    misc/misc_time.lo \
    misc/alloc.lo

## proc2mod:

ngproc2mod_SOURCES = ngproc2mod.c

ngproc2mod_LDADD =  \
	frontend/libfte.la \
	frontend/parser/libparser.la \
	spicelib/parser/libinp.la \
	misc/libmisc.la


## multidec:

ngmultidec_SOURCES = ngmultidec.c

ngmultidec_LDADD =  \
	frontend/libfte.la \
	frontend/parser/libparser.la \
	maths/sparse/libsparse.la \
	misc/libmisc.la

if KLU_WANTED
ngmultidec_LDADD += maths/KLU/libKLU.la
endif


## ngmakeidx:

ngmakeidx_SOURCES = makeidx.c

ngmakeidx_LDADD =  \
	frontend/libfte.la \
	frontend/parser/libparser.la \
	misc/libmisc.la

## create index for online help:

ngspice.idx: ngmakeidx$(EXEEXT) $(srcdir)/ngspice.txt
	./ngmakeidx$(EXEEXT) -o ngspice.idx $(srcdir)/ngspice.txt

endif !SHWIN
endif !WINGUI
endif OLDAPPS

if RELPATH
spinitpath=../lib/ngspice
endif RELPATH
if !RELPATH
spinitpath=$(pkglibdir)
endif !RELPATH

## create spinit at compile time, see
##   (info "(autoconf) Installation Directory Variables")

edit = sed \
	-e 's|@XSPICEINIT[@]|$(XSPICEINIT)|g' \
	-e 's|@pkglibdir[@]|$(spinitpath)|g'

spinit tclspinit : Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	  test -f ./$@.in || srcdir=$(srcdir)/; \
	  $(edit) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

spinit: $(srcdir)/spinit.in
tclspinit: $(srcdir)/tclspinit.in


## General Includes and libraries:

AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_srcdir)/src/include -I$(top_srcdir)/src/spicelib/devices @X_CFLAGS@

AM_CFLAGS = -static

LIBS = @LIBS@ @X_LIBS@ @X_PRE_LIBS@ @X_EXTRA_LIBS@

if WINGUI
LIBS += -lpsapi -lshlwapi
endif

if WINCONSOLE
LIBS += -lpsapi -lshlwapi
endif

if SHWIN
LIBS += -lpsapi -lshlwapi
endif

CLEANFILES = ngspice.idx spinit tclspinit pkgIndex.tcl

MAINTAINERCLEANFILES = Makefile.in

## Tcl:

if TCL_MODULE
lib_LTLIBRARIES = libspice.la
initdata_DATA += tclspinit
pkgIndex_DATA = pkgIndex.tcl
pkgIndexdir = $(pkglibdir)

if TCLWIN
LIBS += -lpsapi -lshlwapi  -ltcl84 -lBLT24
endif

libspice_la_SOURCES = \
	main.c \
	conf.c \
	conf.h \
	ngspice.c \
	tclspice.c

libspice_la_CPPFLAGS = $(AM_CPPFLAGS)  -DSIMULATOR \
	-DTCLSPICE_version="\"$(TCLSPICE_VERSION)\""

libspice_la_LIBADD = \
	@X_LIBS@ \
	frontend/libfte.la

if WINGUI
libspice_la_LIBADD += \
	frontend/wdisp/libwindisp.la
endif

libspice_la_LIBADD += \
	frontend/plotting/libplotting.la

if XSPICE_WANTED
libspice_la_LIBADD += \
	xspice/cm/libcmxsp.la \
	xspice/mif/libmifxsp.la
endif

libspice_la_LIBADD += \
	$(DYNAMIC_DEVICELIBS)

if NDEV_WANTED
libspice_la_LIBADD += \
	spicelib/devices/ndev/libndev.la
endif

if NUMDEV_WANTED
libspice_la_LIBADD += \
	spicelib/devices/nbjt/libnbjt.la   \
	spicelib/devices/nbjt2/libnbjt2.la \
	spicelib/devices/numd/libnumd.la   \
	spicelib/devices/numd2/libnumd2.la \
	spicelib/devices/numos/libnumos.la
endif

libspice_la_LIBADD += \
	spicelib/analysis/libckt.la \
	spicelib/devices/libdev.la

if XSPICE_WANTED
libspice_la_LIBADD += \
	xspice/evt/libevtxsp.la \
	xspice/enh/libenhxsp.la \
	xspice/ipc/libipcxsp.la \
	xspice/idn/libidnxsp.la 
endif
libspice_la_LIBADD += $(XSPICEDLLIBS)

libspice_la_LIBADD += \
	frontend/parser/libparser.la \
	frontend/numparam/libnumparam.la \
	frontend/trannoise/libtrannoise.la \
	spicelib/parser/libinp.la

if CIDER_WANTED
libspice_la_LIBADD += \
	ciderlib/twod/libcidertwod.la     \
	ciderlib/oned/libcideroned.la     \
	ciderlib/input/libciderinput.la   \
	ciderlib/support/libcidersuprt.la
endif

if OSDI_WANTED
libspice_la_LIBADD += osdi/libosdi.la
endif

if KLU_WANTED
libspice_la_LIBADD += maths/KLU/libKLU.la
endif

libspice_la_LIBADD += \
	maths/deriv/libderiv.la \
	maths/cmaths/libcmaths.la \
	maths/misc/libmathmisc.la \
	maths/fft/libmathfft.la \
	maths/poly/libpoly.la \
	maths/ni/libni.la \
	maths/sparse/libsparse.la \
	misc/libmisc.la \
	$(LIBS)

if SP_WANTED
libspice_la_LIBADD += \
    maths/dense/libdense.la
endif

libspice_la_CFLAGS = -shared
libspice_la_LDFLAGS =  -shared -Wl,--version-script=$(srcdir)/tclspice.map
if TCLWIN
libspice_la_LDFLAGS +=  -no-undefined
libspice_la_LDFLAGS += --def spice.def @TCL_LIB_SPEC@
endif
if TCLCYG
libspice_la_LDFLAGS +=  -no-undefined
libspice_la_LDFLAGS += --def spice.def @TCL_LIB_SPEC@
endif

TCLSPICE_VERSION = @VERSION@

pkgIndex_edit = sed \
	-e 's|%LIB_DIR%|$(libdir)|g' \
	-e 's|%VERSION%|$(TCLSPICE_VERSION)|g'

pkgIndex.tcl : Makefile
	rm -f $@ $@.tmp
	srcdir=''; \
	  test -f ./$@.in || srcdir=$(srcdir)/; \
	  $(pkgIndex_edit) $${srcdir}$@.in >$@.tmp
	mv $@.tmp $@

pkgIndex.tcl: $(srcdir)/pkgIndex.tcl.in

endif TCL_MODULE

EXTRA_DIST += spice.def



## ngspice as shared module:

if SHARED_MODULE
lib_LTLIBRARIES = libngspice.la

libngspice_la_SOURCES = \
	conf.c \
	conf.h \
	ngspice.c \
	sharedspice.c

libngspice_la_CPPFLAGS = $(AM_CPPFLAGS)  -DSIMULATOR \
	-DSHAREDSPICE_version="\"$(SHAREDSPICE_VERSION)\""

libngspice_la_LIBADD = \
	@X_LIBS@ \
	frontend/libfte.la

libngspice_la_LIBADD += \
	frontend/plotting/libplotting.la

if XSPICE_WANTED
libngspice_la_LIBADD += \
	xspice/cm/libcmxsp.la \
	xspice/mif/libmifxsp.la
endif

libngspice_la_LIBADD += \
	$(DYNAMIC_DEVICELIBS)

if NDEV_WANTED
libngspice_la_LIBADD += \
	spicelib/devices/ndev/libndev.la
endif

if NUMDEV_WANTED
libngspice_la_LIBADD += \
	spicelib/devices/nbjt/libnbjt.la   \
	spicelib/devices/nbjt2/libnbjt2.la \
	spicelib/devices/numd/libnumd.la   \
	spicelib/devices/numd2/libnumd2.la \
	spicelib/devices/numos/libnumos.la
endif

libngspice_la_LIBADD += \
	spicelib/analysis/libckt.la \
	spicelib/devices/libdev.la

if XSPICE_WANTED
libngspice_la_LIBADD += \
	xspice/evt/libevtxsp.la \
	xspice/enh/libenhxsp.la \
	xspice/ipc/libipcxsp.la \
	xspice/idn/libidnxsp.la 
endif
libngspice_la_LIBADD += $(XSPICEDLLIBS)

libngspice_la_LIBADD += \
	frontend/parser/libparser.la \
	frontend/numparam/libnumparam.la \
	frontend/trannoise/libtrannoise.la \
	spicelib/parser/libinp.la

if CIDER_WANTED
libngspice_la_LIBADD += \
	ciderlib/twod/libcidertwod.la     \
	ciderlib/oned/libcideroned.la     \
	ciderlib/input/libciderinput.la   \
	ciderlib/support/libcidersuprt.la
endif

if OSDI_WANTED
libngspice_la_LIBADD += osdi/libosdi.la
endif

if KLU_WANTED
libngspice_la_LIBADD += maths/KLU/libKLU.la
endif

libngspice_la_LIBADD += \
	maths/deriv/libderiv.la \
	maths/cmaths/libcmaths.la \
	maths/misc/libmathmisc.la \
	maths/fft/libmathfft.la \
	maths/poly/libpoly.la \
	maths/ni/libni.la \
	maths/sparse/libsparse.la \
	misc/libmisc.la \
	$(LIBS) \
	$(OPENMP_CFLAGS)

if SP_WANTED
libngspice_la_LIBADD += \
    maths/dense/libdense.la
endif

libngspice_la_CFLAGS = -shared

libngspice_la_LDFLAGS =  -shared
libngspice_la_LDFLAGS +=  -lstdc++

if SHWIN
libngspice_la_LDFLAGS += -Wl,--output-def=ngspice.def  -Wl,--out-implib=ngspice.dll.a
libngspice_la_LDFLAGS += -no-undefined
EXTRA_DIST += ngspice.def
endif

if SHCYG
libngspice_la_LDFLAGS += -Wl,--output-def=ngspice.def  -Wl,--out-implib=ngspice.dll.a
libngspice_la_LDFLAGS += -no-undefined
EXTRA_DIST += ngspice.def
endif

## if !SHWIN
## if !SHCYG
libngspice_la_LDFLAGS += -version-info $(LT_NGSPICE_CURRENT):$(LT_NGSPICE_REVISION):$(LT_NGSPICE_AGE)
## endif
## endif

SHAREDSPICE_VERSION = @VERSION@

endif SHARED_MODULE
