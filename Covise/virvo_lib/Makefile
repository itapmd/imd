
# Makefile for application library

# All common UNIX compiler flags and options are now in
# this central makefile.

include $(COVISEDIR)/src/make.env.$(ARCHSUFFIX)
include $(COVISEDIR)/src/make.rules

PROJECT = General
APPLICATION = $(VIRVO_LIB_FILE)


OBJECTS = $(OBJDIR)/vvdebugmsg.o  $(OBJDIR)/vvgltools.o  $(OBJDIR)/vvrenderer.o  \
        $(OBJDIR)/vvsoftper.o    $(OBJDIR)/vvtexrend.o    $(OBJDIR)/vvvecmath.o \
        $(OBJDIR)/vvdicom.o     $(OBJDIR)/vvimage.o     $(OBJDIR)/vvsocket.o    \
        $(OBJDIR)/vvsoftvr.o     $(OBJDIR)/vvtokenizer.o  $(OBJDIR)/vvvffile.o \
        $(OBJDIR)/vvdynlib.o   $(OBJDIR)/vvimgrend.o   $(OBJDIR)/vvsocketio.o  \
        $(OBJDIR)/vvsphere.o     $(OBJDIR)/vvtoolshed.o   $(OBJDIR)/vvvoldesc.o \
        $(OBJDIR)/vvfileio.o    $(OBJDIR)/vvprintgl.o   $(OBJDIR)/vvsoftimg.o   \
        $(OBJDIR)/vvstopwatch.o  $(OBJDIR)/vvtransfunc.o   $(OBJDIR)/vvxvid.o $(IDENTOBJ)


all:   $(OBJDIR) $(PLATFORM)/bin/$(VIRVO_LIB_FILE)

$(PLATFORM)/bin/$(VIRVO_LIB_FILE): $(OBJECTS)
	$(MKSHARED2) $(VIRVO_LIB_FILE) $(OBJECTS) $(ADDOBJECTS) $(OGLLIB) $(XVIDLIB) $(XVIDOPTIONS)
	chmod 755 $(VIRVO_LIB_FILE)
	cp $(VIRVO_LIB_FILE) $(PLATFORM)/bin
	@echo VirvoLib created


THISDIR = src/virvo_lib
SRCFILES = $(THISDIR)/*.cpp $(THISDIR)/*.h $(THISDIR)/Makefile
DEVFILES = $(THISDIR)/*.h
BINFILES = $(ARCHSUFFIX)/bin/$(VIRVO_LIB_FILE)

$(OBJDIR)/%.o: %.cpp
	$(LIBPRECC) $<
	$(LIBCC)$< $(LIBCPPFLAGS) $(XVIDOPTIONS) $(XVIDINC)
	$(POSTCC)


# DO NOT DELETE

$(OBJDIR)/vvdebugmsg.o: vvdebugmsg.h
$(OBJDIR)/vvdicom.o: vvdebugmsg.h vvdicom.h vvtoolshed.h
$(OBJDIR)/vvdynlib.o: vvdynlib.h
$(OBJDIR)/vvfileio.o: vvfileio.h vvvoldesc.h vvtoolshed.h vvvecmath.h
$(OBJDIR)/vvfileio.o: vvtransfunc.h vvsllist.h vvvffile.h vvdebugmsg.h
$(OBJDIR)/vvfileio.o: vvtokenizer.h vvdicom.h vvarray.h
$(OBJDIR)/vvgltools.o: vvgltools.h
$(OBJDIR)/vvimage.o: vvimage.h vvtoolshed.h vvdebugmsg.h
$(OBJDIR)/vvimgrend.o: vvdebugmsg.h vvimgrend.h vvrenderer.h vvvecmath.h
$(OBJDIR)/vvimgrend.o: vvvoldesc.h vvtoolshed.h vvtransfunc.h vvsllist.h
$(OBJDIR)/vvimgrend.o: vvimage.h vvsocketio.h vvsocket.h vvsoftimg.h
$(OBJDIR)/vvimgrend.o: vvstopwatch.h
$(OBJDIR)/vvprintgl.o: vvtoolshed.h vvprintgl.h vvdebugmsg.h
$(OBJDIR)/vvrenderer.o: vvrenderer.h vvvecmath.h vvvoldesc.h vvtoolshed.h
$(OBJDIR)/vvrenderer.o: vvtransfunc.h vvsllist.h vvdebugmsg.h vvprintgl.h
$(OBJDIR)/vvsocket.o: vvsocket.h
$(OBJDIR)/vvsocketio.o: vvsocketio.h vvsocket.h vvvoldesc.h vvtoolshed.h
$(OBJDIR)/vvsocketio.o: vvvecmath.h vvtransfunc.h vvsllist.h vvimage.h
$(OBJDIR)/vvsocketio.o: vvdebugmsg.h
$(OBJDIR)/vvsoftimg.o: vvdebugmsg.h vvvecmath.h vvsoftimg.h vvrenderer.h
$(OBJDIR)/vvsoftimg.o: vvvoldesc.h vvtoolshed.h vvtransfunc.h vvsllist.h
$(OBJDIR)/vvsoftper.o: vvdebugmsg.h vvvecmath.h vvsoftper.h vvrenderer.h
$(OBJDIR)/vvsoftper.o: vvvoldesc.h vvtoolshed.h vvtransfunc.h vvsllist.h
$(OBJDIR)/vvsoftper.o: vvsoftimg.h vvsoftvr.h vvimage.h vvstopwatch.h
$(OBJDIR)/vvsoftvr.o: vvdebugmsg.h vvvecmath.h vvsoftimg.h vvrenderer.h
$(OBJDIR)/vvsoftvr.o: vvvoldesc.h vvtoolshed.h vvtransfunc.h vvsllist.h
$(OBJDIR)/vvsoftvr.o: vvsoftvr.h vvimage.h vvstopwatch.h
$(OBJDIR)/vvsphere.o: vvsphere.h vvvecmath.h vvdebugmsg.h vvtoolshed.h
$(OBJDIR)/vvstopwatch.o: vvstopwatch.h
$(OBJDIR)/vvtexrend.o: vvvecmath.h vvdebugmsg.h vvtoolshed.h vvgltools.h
$(OBJDIR)/vvtexrend.o: vvsphere.h vvtexrend.h vvvoldesc.h vvtransfunc.h
$(OBJDIR)/vvtexrend.o: vvsllist.h vvrenderer.h vvstopwatch.h
$(OBJDIR)/vvtokenizer.o: vvtokenizer.h
$(OBJDIR)/vvtoolshed.o: vvtoolshed.h
$(OBJDIR)/vvtransfunc.o: vvdebugmsg.h vvvecmath.h vvtransfunc.h vvtoolshed.h
$(OBJDIR)/vvvecmath.o: vvvecmath.h
$(OBJDIR)/vvvffile.o: vvdebugmsg.h vvvffile.h
$(OBJDIR)/vvvoldesc.o: vvdebugmsg.h vvtoolshed.h vvvecmath.h vvvoldesc.h
$(OBJDIR)/vvvoldesc.o: vvtransfunc.h vvsllist.h
