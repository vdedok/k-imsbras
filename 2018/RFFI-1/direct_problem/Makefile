CFLAGS	= -Wall -g3
LFLAGS	= 
MAKE	= gmake
CC	= g++
SOURCES	= $(wildcard *.cpp)
HEADERS	= $(wildcard *.hpp)
OBJECTS	= $(patsubst %.cpp, %.o, $(SOURCES))
TARGET	= eikonal

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

%.d: %.cpp
	$(SHELL) -ec '$(CC) -MM $(CFLAGS) $< > $@'

all: $(OBJECTS)
	$(CC) -o $(TARGET) $(LFLAGS) $(OBJECTS)

clean:
	find -name "*.o" -o -name "*.d" |xargs rm -rf
	rm -rf $(TARGET)

etags:
	find -name "*.cpp" -o -name "*.c" -o -name "*.hpp" -o -name "*.h" |xargs etags

include $(OBJECTS:.o=.d)

