
PLUGIN_EXT	?= .so
PLUGIN		?= match-vamp-plugin$(PLUGIN_EXT)
CXX		?= g++
CC		?= gcc

HEADERS		:= Finder.h Matcher.h MatchFeeder.h MatchVampPlugin.h Path.h
SOURCES		:= Finder.cpp Matcher.cpp MatchFeeder.cpp MatchVampPlugin.cpp Path.cpp

OBJECTS		:= $(SOURCES:.cpp=.o)

$(PLUGIN):	$(OBJECTS)
		$(CXX) -o $@ $^ $(LDFLAGS)

clean:		
		rm -f $(OBJECTS)

distclean:	clean
		rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)

# DO NOT DELETE

Finder.o: Finder.h Matcher.h
Matcher.o: Matcher.h
MatchFeeder.o: MatchFeeder.h Matcher.h Finder.h
MatchVampPlugin.o: MatchVampPlugin.h Matcher.h MatchFeeder.h Finder.h Path.h
Path.o: Path.h
Finder.o: Matcher.h
MatchFeeder.o: Matcher.h Finder.h
