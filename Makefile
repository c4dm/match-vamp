
CXXFLAGS	:= -I../vamp-plugin-sdk -O3 -Wall
#CXXFLAGS	:= -I../vamp-plugin-sdk -g -Wall -march=pentium4 -msse -msse2 -ffast-math
#CXXFLAGS	:= -I../vamp-plugin-sdk -O3 -Wall -march=pentium4 -msse -msse2 -fomit-frame-pointer -ffast-math

match-vamp-plugin.so:	Finder.o Matcher.o MatchFeeder.o MatchVampPlugin.o Path.o
	g++ -shared $^ -o $@ -L../vamp-plugin-sdk/vamp-sdk -Wl,-Bstatic -lvamp-sdk -Wl,-Bdynamic -lpthread

clean:	
	rm *.o

# DO NOT DELETE

Finder.o: Finder.h Matcher.h
Matcher.o: Matcher.h Finder.h
MatchFeeder.o: MatchFeeder.h Matcher.h Finder.h
MatchVampPlugin.o: MatchVampPlugin.h Matcher.h MatchFeeder.h Finder.h Path.h
Path.o: Path.h
