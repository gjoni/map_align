
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

all: $(OBJDIR) map_align mp_super

$(OBJDIR):
	mkdir -p $(OBJDIR)

map_align: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o map_align $(filter-out obj/mp_super.o, $(OBJS)) $(CXXLIBS)

mp_super: $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCDIRS) $(LIBDIRS) -o mp_super obj/mp_super.o obj/Chain.o obj/Residue.o obj/Atom.o obj/kdtree.o obj/TMalign.o obj/Kabsch.o $(CXXLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCDIRS) -c -o $@ $<

clean:
	rm $(OBJS) map_align mp_super
	rmdir $(OBJDIR)
