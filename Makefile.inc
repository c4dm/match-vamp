
PLUGIN_EXT	?= .so
PLUGIN		?= match-vamp-plugin$(PLUGIN_EXT)
CXX		?= g++
CC		?= gcc

HEADERS		:= DistanceMetric.h Finder.h Matcher.h MatchFeeder.h MatchFeatureFeeder.h MatchVampPlugin.h Path.h
SOURCES		:= DistanceMetric.cpp Finder.cpp Matcher.cpp MatchFeeder.cpp MatchFeatureFeeder.cpp MatchVampPlugin.cpp Path.cpp

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

DistanceMetric.o: DistanceMetric.h
Finder.o: Finder.h Matcher.h DistanceMetric.h
Matcher.o: Matcher.h DistanceMetric.h
MatchFeeder.o: MatchFeeder.h Matcher.h DistanceMetric.h Finder.h
MatchFeatureFeeder.o: MatchFeatureFeeder.h Matcher.h DistanceMetric.h
MatchFeatureFeeder.o: Finder.h
MatchVampPlugin.o: MatchVampPlugin.h Matcher.h DistanceMetric.h MatchFeeder.h
MatchVampPlugin.o: Finder.h Path.h
Path.o: Path.h
Finder.o: Matcher.h DistanceMetric.h
Matcher.o: DistanceMetric.h
MatchFeeder.o: Matcher.h DistanceMetric.h Finder.h
MatchFeatureFeeder.o: Matcher.h DistanceMetric.h Finder.h
MatchVampPlugin.o: Matcher.h DistanceMetric.h
