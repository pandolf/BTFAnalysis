CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


<<<<<<< HEAD
all: positionAnalysis drawPositionAnalysis hodoscopeAnalysis
=======
all: positionAnalysis drawPositionAnalysis calibrateCeF3
>>>>>>> 22c506bb8f0fcd0ad32a32a76c963fdbdbf17226



positionAnalysis: positionAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o positionAnalysis positionAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

drawPositionAnalysis: drawPositionAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o drawPositionAnalysis drawPositionAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)

calibrateCeF3: calibrateCeF3.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateCeF3 calibrateCeF3.cpp $(ROOTFLAG) $(EXTRALIBS)

hodoscopeAnalysis: hodoscopeAnalysis.cpp
	$(CC) -Wall $(INCLUDES) -o hodoscopeAnalysis hodoscopeAnalysis.cpp $(ROOTFLAG) $(EXTRALIBS)



clean:
	rm positionAnalysis
	rm drawPositionAnalysis 
	rm hodoscopeAnalysis
