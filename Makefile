
CXXFLAGS	:= -fPIC -ffast-math -O3 -Wall

match-vamp-plugin.so:	Finder.o Matcher.o MatchFeeder.o MatchVampPlugin.o Path.o
	g++ -shared $^ -o $@ -Wl,-Bstatic -lvamp-sdk -Wl,-Bdynamic -lpthread -Wl,--version-script=vamp-plugin.map

clean:	
	rm *.o

# DO NOT DELETE

Finder.o: Finder.h Matcher.h
Matcher.o: Matcher.h Finder.h
MatchFeeder.o: MatchFeeder.h Matcher.h Finder.h
MatchVampPlugin.o: MatchVampPlugin.h Matcher.h MatchFeeder.h Finder.h Path.h
Path.o: Path.h
