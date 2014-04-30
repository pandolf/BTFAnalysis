CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


all: positionAnalysis drawPositionAnalysis

positionAnalysis: positionAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

drawPositionAnalysis: drawPositionAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o drawPositionAnalysis drawPositionAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)




clean:
#	rm *.o
	rm positionAnalysis
	rm drawPositionAnalysis 
