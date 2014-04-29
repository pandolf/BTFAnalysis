CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib  -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit




calibrateFibresWithCosmics: calibrateFibresWithCosmics.cpp
	$(CC) -Wall $(INCLUDES) -o calibrateFibresWithCosmics calibrateFibresWithCosmics.cpp $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
