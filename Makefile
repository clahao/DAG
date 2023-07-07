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
TARGETS_SSSP_SINGLE  =  sssp_single
TARGETS_SSSP_SCHEDULE  =  sssp_schedule
TARGETS_SSSP_PARALLEL_SYNC  =  sssp_parallel_sync
TARGETS_SSSP_PARALLEL_ASYNC  =  sssp_parallel_async
TARGETS_SSSP_PARALLEL_SCC  =  sssp_parallel_scc
TARGETS_SSWP_SCC   = sswp_scc
TARGETS_SSWP_ASYNC  =  sswp_async
TARGETS_SSWP_SYNC  =  sswp_sync
TARGETS_SSWP_SINGLE  =  sswp_single
TARGETS_CC_SCC   = cc_scc
TARGETS_CC_ASYNC  =  cc_async
TARGETS_CC_SYNC  =  cc_sync
TARGETS_CC_SINGLE  =  cc_single
TARGETS_PR_SCC   = pr_scc
TARGETS_PR_ASYNC  =  pr_async
TARGETS_PR_SYNC  =  pr_sync
TARGETS_PR_SINGLE  =  pr_single
TARGETS_DENSITY  =  community_density

all : $(TARGETS_SSSP_SCC) $(TARGETS_SSSP_ASYNC) $(TARGETS_SSSP_SYNC) $(TARGETS_SSSP_SINGLE) $(TARGETS_SSSP_SCHEDULE) $(TARGETS_SSSP_PARALLEL_SYNC) $(TARGETS_SSSP_PARALLEL_SCC)\
	  $(TARGETS_SSWP_SCC) $(TARGETS_SSWP_ASYNC) $(TARGETS_SSWP_SYNC) $(TARGETS_SSWP_SINGLE)\
	  $(TARGETS_CC_SCC) $(TARGETS_CC_ASYNC) $(TARGETS_CC_SYNC) $(TARGETS_CC_SINGLE)\
	  $(TARGETS_PR_SCC) $(TARGETS_PR_ASYNC) $(TARGETS_PR_SYNC) $(TARGETS_PR_SINGLE)\
	  $(TARGETS_DENSITY)

sssp : $(TARGETS_SSSP_SCC) $(TARGETS_SSSP_ASYNC) $(TARGETS_SSSP_SYNC)  $(TARGETS_SSSP_SCHEDULE)
sswp : $(TARGETS_SSWP_SCC) $(TARGETS_SSWP_ASYNC) $(TARGETS_SSWP_SYNC)
cc : $(TARGETS_CC_SCC) $(TARGETS_CC_ASYNC) $(TARGETS_CC_SYNC)
pr : $(TARGETS_PR_SCC) $(TARGETS_PR_ASYNC) $(TARGETS_PR_SYNC)
single: $(TARGETS_SSSP_SINGLE) $(TARGETS_SSWP_SINGLE) $(TARGETS_CC_SINGLE) $(TARGETS_PR_SINGLE)
parallel: $(TARGETS_SSSP_PARALLEL_SYNC) $(TARGETS_SSSP_PARALLEL_ASYNC) $(TARGETS_SSSP_PARALLEL_SCC)
test: $(TARGETS_DENSITY)

#@g++ -o $@ $< -w
$(TARGETS_SSSP_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSSP_ASYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSSP_SYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSSP_SINGLE): %: %.cpp
	@g++ -mavx512f -o $@ $< -w

$(TARGETS_SSSP_SCHEDULE): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d
#@g++ -mavx512f -o $@ $< -w
$(TARGETS_SSSP_PARALLEL_SYNC): %: %.cpp
	@g++ -mavx512f -o $@ $< -w

$(TARGETS_SSSP_PARALLEL_ASYNC): %: %.cpp
	@g++ -mavx512f -o $@ $< -w

$(TARGETS_SSSP_PARALLEL_SCC): %: %.cpp
	@g++ -mavx512f -o $@ $< -w

$(TARGETS_SSWP_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	      -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSWP_ASYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSWP_SYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_SSWP_SINGLE): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_CC_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_CC_ASYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_CC_SYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_CC_SINGLE): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_PR_SCC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_PR_ASYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_PR_SYNC): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_PR_SINGLE): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

$(TARGETS_DENSITY): %: %.cpp
	@$(LINK.cc) -MD -o $@ $< $(LDLIBS) $(CFLAGS) $(LDFLAGS) -w
	@cp $*.d .$*.P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> .$*.P; \
	  rm -f $*.d

.PHONY: clean
clean:
	@$(RM) $(TARGETS_SSSP_SCC) $(TARGETS_SSSP_SCC:%=.%.P) $(TARGETS_SSSP_ASYNC) $(TARGETS_SSSP_ASYNC:%=.%.P) $(TARGETS_SSSP_SYNC) $(TARGETS_SSSP_SYNC:%=.%.P) \
	 $(TARGETS_SSSP_SINGLE) $(TARGETS_SSSP_SINGLE:%=.%.P) $(TARGETS_SSSP_SCHEDULE) $(TARGETS_SSSP_SCHEDULE:%=.%.P) \
	 $(TARGETS_SSSP_PARALLEL_SYNC) $(TARGETS_SSSP_PARALLEL_SYNC:%=.%.P) $(TARGETS_SSSP_PARALLEL_SCC) $(TARGETS_SSSP_PARALLEL_SCC:%=.%.P)\
     $(TARGETS_SSWP_SCC) $(TARGETS_SSWP_SCC:%=.%.P) $(TARGETS_SSWP_ASYNC) $(TARGETS_SSWP_ASYNC:%=.%.P) $(TARGETS_SSWP_SYNC) $(TARGETS_SSWP_SYNC:%=.%.P) \
	 $(TARGETS_CC_SCC) $(TARGETS_CC_SCC:%=.%.P) $(TARGETS_CC_ASYNC) $(TARGETS_CC_ASYNC:%=.%.P) $(TARGETS_CC_SYNC) $(TARGETS_CC_SYNC:%=.%.P)	\
	 $(TARGETS_PR_SCC) $(TARGETS_PR_SCC:%=.%.P) $(TARGETS_PR_ASYNC) $(TARGETS_PR_ASYNC:%=.%.P) $(TARGETS_PR_SYNC) $(TARGETS_PR_SYNC:%=.%.P)


-include .*.P

