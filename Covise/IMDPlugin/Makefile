###########################################################
#
#  Usage:
#
#    source $COVISEDIR/.covise $COVISEDIR
#    make
#    copy libIMDPlugin.so to $COVISEDIR/$ARCHSUFFIX/bin
#
###########################################################

COVISEDIR  = /common/linux/paket/Covise/covise_5.3.1/covise
ARCHSUFFIX = gcc3

export COVISEDIR ARCHSUFFIX

APPLICATION = IMDPlugin
PROJECT = SFB382

#VIRVODIR = $(COVISEDIR)/src/renderer/COVER/plugins/Volume
VIRVODIR = ../virvo_lib
VIRVOINC = -I$(VIRVODIR)

IMDDIR = ../../src
IMDINC = -I$(IMDDIR)

## plugins need COVER headers to compile
ADD_INC =  $(ARTOOLKITINC) $(ARTOOLKITOPTIONS) $(COVISE_BRANCH) $(CXX_STANDARD) -I$(COVISEDIR)/src -I$(COVISEDIR)/src/renderer/COVER -I$(COVISEDIR)/src/renderer/COVER/vrui -I$(COVISEDIR)/src/covise_util -I$(PERFORMERHOME)/include -I$(OGLHOME)/include $(XINCLUDE) $(VIRVOINC) $(IMDINC)

## and they build a special object
#ADD_CLEAN = $(PLATFORM)/bin/lib$(APPLICATION).so
ADD_CLEAN = lib$(APPLICATION).so

## add other objects here if you have more than
OBJECTS = $(OBJDIR)/$(APPLICATION).o $(OBJDIR)/AtomScenegraph.o \
  $(OBJDIR)/LimitsDialog.o

# CATEGORY: no space after '=' !!! Plugins have no category !!!
CATEGORY =

# All common UNIX compiler flags and options are in
# this central makefile.
include $(COVISEDIR)/src/make.env.$(ARCHSUFFIX)

## our main target
#plugin: $(PLATFORM)/bin/lib$(APPLICATION).so
plugin: lib$(APPLICATION).so

# All common rules are in this central makefile.
include $(COVISEDIR)/src/make.rules

# -------------------------------------------------------------------------

#$(PLATFORM)/bin/lib$(APPLICATION).so: $(OBJECTS)
#	$(MKSHARED) $(PLATFORM)/bin/lib$(APPLICATION).so $(OBJECTS) $(VIRVOLIB)
lib$(APPLICATION).so: $(OBJECTS)
	$(MKSHARED) lib$(APPLICATION).so $(OBJECTS) $(VIRVOLIB)

$(OBJDIR)/$(APPLICATION).o: $(OBJDIR) $(APPLICATION).cpp $(APPLICATION).h
	$(CC) $(OPTIONS) $(OPT) $(COVISEINC) $(ADD_INC) -c -o $(OBJDIR)/$(APPLICATION).o $(APPLICATION).cpp

$(OBJDIR)/AtomScenegraph.o: AtomScenegraph.cpp
	$(CC) $(CPPFLAGS) -c AtomScenegraph.cpp -o $(OBJDIR)/AtomScenegraph.o

$(OBJDIR)/LimitsDialog.o: LimitsDialog.cpp
	$(CC) $(CPPFLAGS) -c LimitsDialog.cpp -o $(OBJDIR)/LimitsDialog.o

$(OBJDIR)/vvsocket.o: $(VIRVODIR)/vvsocket.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvsocket.cpp -o $(OBJDIR)/vvsocket.o

$(OBJDIR)/vvsocketio.o: $(VIRVODIR)/vvsocketio.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvsocketio.cpp -o $(OBJDIR)/vvsocketio.o

$(OBJDIR)/vvtoolshed.o: $(VIRVODIR)/vvtoolshed.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvtoolshed.cpp -o $(OBJDIR)/vvtoolshed.o

$(OBJDIR)/vvvoldesc.o: $(VIRVODIR)/vvvoldesc.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvvoldesc.cpp -o $(OBJDIR)/vvvoldesc.o

$(OBJDIR)/vvimage.o: $(VIRVODIR)/vvimage.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvimage.cpp -o $(OBJDIR)/vvimage.o

$(OBJDIR)/vvvecmath.o: $(VIRVODIR)/vvvecmath.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvvecmath.cpp -o $(OBJDIR)/vvvecmath.o

$(OBJDIR)/vvdebugmsg.o: $(VIRVODIR)/vvdebugmsg.cpp
	$(CC) $(CPPFLAGS) -c $(VIRVODIR)/vvdebugmsg.cpp -o $(OBJDIR)/vvdebugmsg.o
# DO NOT DELETE

$(OBJDIR)/AtomScenegraph.o: AtomScenegraph.h
$(OBJDIR)/IMDPlugin.o: IMDPlugin.h AtomScenegraph.h LimitsDialog.h
$(OBJDIR)/LimitsDialog.o: LimitsDialog.h IMDPlugin.h AtomScenegraph.h
