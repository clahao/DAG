CXXFLAGS += -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
	    -Wformat=2 -Winit-self -Wmissing-declarations \
	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
	    -Wswitch-default -Wuninitialized -Wwrite-strings \
	    -Wno-sign-conversion -Wno-unused-function \
            -Wno-missing-declarations \
            -fopenmp -std=c++14 -mcx16 -O3 -DNDEBUG
LDLIBS   += -ltcmalloc_minimal -lnuma
CFLAGS += -I /home/lyj/package/boost_1_58_0/include
LDFLAGS += -L /home/lyj/package/boost_1_58_0/lib
TARGETS   = sssp_scc
TARGETS1  =  sssp_async
TARGETS2  =  sssp_sync
all : $(TARGETS) $(TARGETS1) $(TARGETS2)

$(TARGETS): %: %.cpp
	$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS1): %: %.cpp
	$(LINK.cc) -MD -o $@ $< -w

$(TARGETS2): %: %.cpp
	$(LINK.cc) -MD -o $@ $< -w

.PHONY: clean
clean:
	$(RM) $(TARGETS) $(TARGETS:%=.%.P)

-include .*.P

