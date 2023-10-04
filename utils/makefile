CXX=clang++-12
CXXFLAGS= -I .. -D DSA_TARGET_LINUX -std=c++20 -stdlib=libc++ -pthread -mavx2

# Project files
SRCDIR = .
SHSRCS = aa.cc cdn.cc dna.cc polymer.cc
SRCS = $(SHSRCS) utils.cc
OBJS = $(SRCS:.cc=.o)
DEPS = $(SRCS:.cc=.d)
EXE = dsa-util

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

$(DBGDIR)/%.o: $(SRCDIR)/%.cc makefile
	$(CXX) $(CXXFLAGS) $(DBGCXXFLAGS) -MMD -MP -c $< -o $@

# Release build
release: $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CXX) $(CXXFLAGS) $(RELCXXFLAGS) $^ -o $(RELEXE)

-include $(RELDEPS)

$(RELDIR)/%.o: $(SRCDIR)/%.cc makefile
	$(CXX) $(CXXFLAGS) $(RELCXXFLAGS) -MMD -MP -c $< -o $@

# Misc rules
doc:
	@mkdir -p doc
	doxygen Doxyfile

prep:
	@for src in $(SHSRCS) ; do if [ ! -f $${src} ] ; then ln -s ../$${src} ./$${src} ; fi ; done 
	@mkdir -p $(DBGDIR) $(RELDIR)

remake: clean all

clean:
	rm -rf doc/html
	rm -f $(DBGEXE) $(DBGOBJS) $(DBGDEPS) $(RELEXE) $(RELOBJS) $(RELDEPS)
