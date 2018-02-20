
CXX=g++
CXXFLAGS += -Wall -Wno-unused-result -pedantic -Ofast -march=native -std=c++0x -g -ggdb -fopenmp
CXXLIBS =
INCDIRS = -I./src
LIBDIRS =

OBJDIR = obj
SRCDIR = src

DEPS = $(shell find $(SRCDIR) -name '*.h')
SRCS = $(shell find $(SRCDIR) -name '*.cpp')
OBJS = $(patsubst $(SRCDIR)%.cpp, $(OBJDIR)%.o, $(SRCS))

all: $(OBJDIR) map_align

$(OBJDIR):
	mkdir -p $(OBJDIR)

map_align: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o map_align $(OBJS) $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) map_align
	rmdir $(OBJDIR)
