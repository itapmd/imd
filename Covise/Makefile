##########################################################
## Description: Makefile Template for COVISE modules    ##
## (C) 1995 Computer Center University of Stuttgart     ##
##            Allmandring 30  70550 Stuttgart           ##  
## Author: Uwe Woessner, Andreas Wierse  Date: 29.08.97 ##
##########################################################

ARCH=${HOSTTYPE}
ifneq (,$(strip $(findstring linux,${HOSTTYPE})))
ARCH = linux
CC = /opt/experimental/bin/gcc
LIBDIR = -L/opt/experimental/lib
endif
ifneq (,$(strip $(findstring iris,${HOSTTYPE})))
ARCH = irix
endif
COVISEDIR = /common/$(ARCH)/paket/Covise/covise_5.2.3_HLRS_SNAP_0309/covise
OBJDIR = .

# Name of the module executable:
APPLICATION = ReadIMD

# Name of the module category:
CATEGORY = IO_Module

# Needed for creating distribution tar-files:
PROJECT	= sfb382

VIRVODIR = $(COVISEDIR)/src/renderer/COVER/plugins/Volume
ADD_INC = -I$(VIRVODIR) \
          -I$(COVISEDIR)/src/covise_api \
	  -I$(COVISEDIR)/src/appl_lib \
	  -I$(COVISEDIR)/src/covise_lib
COVISE_LIBDIR = -L$(COVISEDIR)/$(ARCH)/bin  
COVISEAPI = -lcovise_api
APPLLIB = -lappl
COVISELIB = -lcovise
OBJECTS = $(OBJDIR)/vvtokenizer.o $(OBJDIR)/vvtoolshed.o
          
# $(APPLICATION) is also the name of the main source and include file!!!

##########################################################

# modify this for your own path if necessary
# COVISEDIR = ../../../..

# All common UNIX compiler flags and options are in
# this central makefile.

# include $(COVISEDIR)/src/make.env.$(ARCHSUFFIX)
#include $(COVISEDIR)/src/make.env

# if you want to modify settings of variables that are set in the
# central makefile, do this here:


# All common rules are in this central makefile.

# include $(COVISEDIR)/src/make.rules

############### link rule for application module #########
# can be modified (e.g. for dependencies or libraries)
##########################################################

#$(PLATFORM)/bin/$(CATEGORY)/$(APPLICATION): $(OBJDIR)/$(APPLICATION).o $(OBJECTS)
$(APPLICATION): $(OBJDIR)/$(APPLICATION).o $(OBJECTS)
	  $(CC) $(CPPFLAGS) -o $@ $(OBJECTS) $(OBJDIR)/$(APPLICATION).o \
	  $(COVISE_LIBDIR) $(COVISEAPI) $(APPLLIB) $(COVISELIB) $(LIBDIR) -lm

############### general object rule ######################
# can be modified (e.g. for dependencies)
##########################################################

$(OBJDIR)/$(APPLICATION).o: $(APPLICATION).cpp $(APPLICATION).h
	$(CC) $(CPPFLAGS) $(ADD_INC) -c $(APPLICATION).cpp -o $(OBJDIR)/$(APPLICATION).o $(LIBDIR)
	$(POSTCC) 

############### module specific rules ####################
# rules for objects in $(OBJECTS)
##########################################################

$(OBJDIR)/vvtokenizer.o: $(VIRVODIR)/vvtokenizer.cpp
	$(CC) $(CPPFLAGS) $(ADD_INC) -c $(VIRVODIR)/vvtokenizer.cpp -o $(OBJDIR)/vvtokenizer.o $(LIBDIR)
	$(POSTCC)

$(OBJDIR)/vvtoolshed.o: $(VIRVODIR)/vvtoolshed.cpp
	$(CC) $(CPPFLAGS) $(ADD_INC) -c $(VIRVODIR)/vvtoolshed.cpp -o $(OBJDIR)/vvtoolshed.o $(LIBDIR)
	$(POSTCC)

##########################################################
# clean: remove .o files
##########################################################

clean:
	rm -f *.o *~ ReadIMD
	rm -rf ii_*

# DO NOT DELETE
