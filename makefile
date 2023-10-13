CXX=clang++
CXXFLAGS= -D DSA_TARGET_LINUX -std=c++20 -mavx2 -stdlib=libc++ -pthread

# Project files
SRCDIR = .
SRCS = aa.cc abs.cc align.cc cdn.cc dna.cc help.cc io.cc main.cc mainfunctions.cc params.cc polymer.cc umi.cc tests.cc
OBJS = $(SRCS:.cc=.o)
DEPS = $(SRCS:.cc=.d)
EXE = dsa

# Debug settings
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGDEPS = $(addprefix $(DBGDIR)/, $(DEPS))
DBGCXXFLAGS = -g -O0 -DDEBUG

# Release settings
RELDIR = bin
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELDEPS = $(addprefix $(RELDIR)/, $(DEPS))
RELCXXFLAGS = -O3 -DNDEBUG

.PHONY: all clean debug doc prep release remake

# Default to release build
all: prep release

# Debug build
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CXX) $(CXXFLAGS) $(DBGCXXFLAGS) $^ -o $(DBGEXE)

-include $(DBGDEPS)

$(DBGDIR)/%.o: $(SRCDIR)/%.cc Makefile
	$(CXX) $(CXXFLAGS) $(DBGCXXFLAGS) -MMD -MP -c $< -o $@

# Release build
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CXX) $(CXXFLAGS) $(RELCXXFLAGS) $^ -o $(RELEXE)

-include $(RELDEPS)

$(RELDIR)/%.o: $(SRCDIR)/%.cc Makefile
	$(CXX) $(CXXFLAGS) $(RELCXXFLAGS) -MMD -MP -c $< -o $@

# Misc rules
doc:
	@mkdir -p doc
	doxygen Doxyfile

prep:
	@mkdir -p $(DBGDIR) $(RELDIR)

remake: clean all

clean:
	rm -rf doc/html
	rm -f $(DBGEXE) $(DBGOBJS) $(DBGDEPS) $(RELEXE) $(RELOBJS) $(RELDEPS)
