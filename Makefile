PROGRAM = gotoh
C_FILES := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp, %.o, $(C_FILES))
CC = g++
CFLAGS = -Wall -pedantic -Ilib -m64 -ffast-math -ftree-vectorize -mssse3 -O3 
LDFLAGS =

all: $(PROGRAM)

$(PROGRAM): .depend $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

depend: .depend

.depend: cmd = gcc -MM -MF depend $(var); cat depend >> .depend;
.depend:
	@echo "Generating dependencies..."
	@$(foreach var, $(C_FILES), $(cmd))
	@rm -f depend

-include .depend

# These are the pattern matching rules. In addition to the automatic
# variables used here, the variable $* that matches whatever % stands for
# can be useful in special cases.
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%: %.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f .depend *.o

.PHONY: clean depend