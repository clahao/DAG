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
TARGETS_SSSP_SCC   = sssp_scc
TARGETS_SSSP_ASYNC  =  sssp_async
TARGETS_SSSP_SYNC  =  sssp_sync
TARGETS_SSWP_SCC   = sswp_scc
TARGETS_SSWP_ASYNC  =  sswp_async
TARGETS_SSWP_SYNC  =  sswp_sync
all : $(TARGETS_SSSP_SCC) $(TARGETS_SSSP_ASYNC) $(TARGETS_SSSP_SYNC) $(TARGETS_SSWP_SCC) $(TARGETS_SSWP_ASYNC) $(TARGETS_SSWP_SYNC)

$(TARGETS_SSSP_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSSP_ASYNC): %: %.cpp
	@g++ -o $@ $< -w

$(TARGETS_SSSP_SYNC): %: %.cpp
	@g++ -o $@ $< -w

$(TARGETS_SSWP_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSWP_ASYNC): %: %.cpp
	@g++ -o $@ $< -w

$(TARGETS_SSWP_SYNC): %: %.cpp
	@g++ -o $@ $< -w

.PHONY: clean
clean:
	@$(RM) $(TARGETS) $(TARGETS:%=.%.P)

-include .*.P

