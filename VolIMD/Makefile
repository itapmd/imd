#
# Use GNU make!!!
#
# Makefile for VolIMD
# author ....................  Roland Niemeier
# version ...................  0.1
# date of last change .......  08/29/96
# description of tool VolIMD: interface for direct volume visualization
# of scalar fields. The application development is for the interactive 
# visualization of volume data sets integrating computer simulations. 
# It is based on Open Inventor, volpack, xforms and a simulation package. 
# Virtual Reality (VR) with direct volume rendering is intended, 
# but actually this is just a work in progress.

# object, source and header files
OBJ = viewer.o render.o sockutil.o client.o menu.o guiforms.o main.o 
SRC = viewer.C render.C menu.C menu.inc main.C ImgRenderAction.C \
	guiforms.C sockutil.C client.c
HDR = viewer.h render.h menu.h imginfo.h global.h ImgRenderAction.h \
	volrend.h guiforms.h tokens.h sockutil.h client.h

# libraries
STDLIBS = -lInventor -lInventorXt -lSgm -lGLU -lGL -lX11 -lXt -lXm -lXext -lm
LIBS    = -L. -lImgViewer -lforms -lvolpack

# compiler and flags
CC       = CC -32 -O            # compiles and runs on juliet and godot
# CC       = CC -n32 -mips3 -O  # compiles and runs on juliet only

.C.o:
	$(CC) -c $<

# all: VolIMD
VolIMD: $(OBJ) libImgViewer.a 
	$(CC) -o VolIMD $(OBJ) $(LIBS) $(STDLIBS) 

# object files
sockutil.o: sockutil.C sockutil.h
	$(CC) -c sockutil.C 
client.o: client.c sockutil.o tokens.h
	$(CC) -c client.c sockutil.o 
viewer.o: viewer.C render.h global.h ImgExaminerViewer.h ImgRenderAction.h 
render.o: render.C viewer.h global.h guiforms.h ImgRenderAction.h 
guiforms.o: guiforms.C global.h render.h
menu.o: menu.C menu.inc global.h render.h ImgRenderAction.h client.h tokens.h
main.o: main.C menu.h viewer.h render.h global.h ImgExaminerViewer.h

# ImgViewerLibrary

libImgViewer.a: libImgViewer.a(InvPixmapButton.o)  \
	libImgViewer.a(ImgViewer.o)                \
	libImgViewer.a(ImgFullViewer.o)            \
	libImgViewer.a(ImgExaminerViewer.o)        \
	libImgViewer.a(ImgRenderAction.o)

libImgViewer.a(ImgRenderAction.o): ImgRenderAction.C ImgRenderAction.h

libImgViewer.a(InvPixmapButton.o):InvPixmapButton.C InvPixmapButton.h 

libImgViewer.a(ImgViewer.o): ImgViewer.C ImgViewer.h ImgRenderAction.h

libImgViewer.a(ImgFullViewer.o): ImgFullViewer.C ImgFullViewer.h \
	ImgViewer.h InvPixmapButton.h ImgRenderAction.h

libImgViewer.a(ImgExaminerViewer.o): ImgExaminerViewer.C \
	ImgExaminerViewer.h ImgFullViewer.h ImgViewer.h ImgRenderAction.h

