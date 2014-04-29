CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit




bgoCalibration: bgoCalibration.cpp
	$(CC) -Wall -Wl,--no-as-needed $(INCLUDES) -o bgoCalibration bgoCalibration.cpp $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
