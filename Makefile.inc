
PLUGIN_EXT	?= .so
PLUGIN		?= match-vamp-plugin$(PLUGIN_EXT)
CXX		?= g++
CC		?= gcc

HEADERS		:= $(wildcard src/*.h)
SOURCES		:= $(wildcard src/*.cpp)

OBJECTS		:= $(SOURCES:.cpp=.o)

$(PLUGIN):	$(OBJECTS)
		$(CXX) -o $@ $^ $(LDFLAGS) $(PLUGIN_LDFLAGS)

clean:		
		rm -f $(OBJECTS)

distclean:	clean
		rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)
# DO NOT DELETE

src/DistanceMetric.o: src/DistanceMetric.h
src/Path.o: src/Path.h
src/MatchFeatureFeeder.o: src/MatchFeatureFeeder.h src/Matcher.h
src/MatchFeatureFeeder.o: src/DistanceMetric.h src/Finder.h
src/FeatureExtractor.o: src/FeatureExtractor.h
src/Finder.o: src/Finder.h src/Matcher.h src/DistanceMetric.h src/Path.h
src/Matcher.o: src/Matcher.h src/DistanceMetric.h
src/MatchVampPlugin.o: src/MatchVampPlugin.h src/Matcher.h
src/MatchVampPlugin.o: src/DistanceMetric.h src/FeatureExtractor.h
src/MatchVampPlugin.o: src/MatchFeatureFeeder.h src/Finder.h src/Path.h
src/MatchFeatureFeeder.o: src/Matcher.h src/DistanceMetric.h src/Finder.h
src/Finder.o: src/Matcher.h src/DistanceMetric.h
src/Matcher.o: src/DistanceMetric.h
src/MatchVampPlugin.o: src/Matcher.h src/DistanceMetric.h
src/MatchVampPlugin.o: src/FeatureExtractor.h
