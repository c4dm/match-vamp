
PLUGIN_EXT	?= .so
PLUGIN		?= match-vamp-plugin$(PLUGIN_EXT)
CXX		?= g++
CC		?= gcc

CXXFLAGS	:= -Isrc -DUSE_COMPACT_TYPES $(CXXFLAGS)

TEST_LDFLAGS	+= -lboost_unit_test_framework

HEADERS		:= $(wildcard src/*.h)
SOURCES		:= $(wildcard src/*.cpp)
OBJECTS		:= $(SOURCES:.cpp=.o)

TEST_SOURCES	:= $(wildcard test/*.cpp)
TEST_OBJECTS	:= $(TEST_SOURCES:.cpp=.o)
TESTS		:= $(TEST_OBJECTS:.o=)

all:		$(PLUGIN) $(TESTS)
		for t in $(TESTS); do echo "Running $$t"; ./"$$t" || exit 1; done

$(PLUGIN):	$(OBJECTS)
		$(CXX) -o $@ $^ $(LDFLAGS) $(PLUGIN_LDFLAGS)

test/%:	test/%.o $(OBJECTS)
		$(CXX) -o $@ $^ $(LDFLAGS) $(TEST_LDFLAGS)

clean:		
		rm -f $(OBJECTS) $(TEST_OBJECTS) $(TESTS)

distclean:	clean
		rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS) $(TEST_SOURCES)
	
# DO NOT DELETE

src/DistanceMetric.o: src/DistanceMetric.h
src/Path.o: src/Path.h
src/FeatureConditioner.o: src/FeatureConditioner.h
src/MatchFeatureFeeder.o: src/MatchFeatureFeeder.h src/Matcher.h
src/MatchFeatureFeeder.o: src/DistanceMetric.h src/MatchTypes.h src/Finder.h
src/FeatureExtractor.o: src/FeatureExtractor.h
src/Finder.o: src/Finder.h src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/Finder.o: src/Path.h
src/Matcher.o: src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/MatchPipeline.o: src/MatchPipeline.h src/Matcher.h src/DistanceMetric.h
src/MatchPipeline.o: src/MatchTypes.h src/Finder.h src/FeatureExtractor.h
src/MatchPipeline.o: src/FeatureConditioner.h src/MatchFeatureFeeder.h
src/MatchVampPlugin.o: src/MatchVampPlugin.h src/MatchPipeline.h
src/MatchVampPlugin.o: src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/MatchVampPlugin.o: src/Finder.h src/FeatureExtractor.h
src/MatchVampPlugin.o: src/FeatureConditioner.h src/MatchFeatureFeeder.h
src/MatchVampPlugin.o: src/Path.h
src/MatchFeatureFeeder.o: src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/MatchFeatureFeeder.o: src/Finder.h
src/Finder.o: src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/Matcher.o: src/DistanceMetric.h src/MatchTypes.h
src/MatchPipeline.o: src/Matcher.h src/DistanceMetric.h src/MatchTypes.h
src/MatchPipeline.o: src/Finder.h src/FeatureExtractor.h
src/MatchPipeline.o: src/FeatureConditioner.h src/MatchFeatureFeeder.h
src/MatchVampPlugin.o: src/MatchPipeline.h src/Matcher.h src/DistanceMetric.h
src/MatchVampPlugin.o: src/MatchTypes.h src/Finder.h src/FeatureExtractor.h
src/MatchVampPlugin.o: src/FeatureConditioner.h src/MatchFeatureFeeder.h
test/TestFeatureConditioner.o: src/FeatureConditioner.h
test/TestDistanceMetric.o: src/DistanceMetric.h
test/TestFeatureExtractor.o: src/FeatureExtractor.h
