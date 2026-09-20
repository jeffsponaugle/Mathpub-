# Build geng variants that generate connected planar graphs by pruning.
#
#   geng_coplanar        : complement trick (min-degree augmentation) -- the fast one
#   geng_coplanar_instr  : same, with geng's INSTRUMENT statistics on stderr
#   geng_direct          : baseline: plain planarity PREPRUNE with geng -c
#
# Requires a configured+built nauty tree in $(NAUTY) (./configure && make).

NAUTY   ?= ../nauty2_9_3
CC      ?= gcc
CFLAGS  ?= -O3 -march=native -DNDEBUG
W1       = -DMAXN=WORDSIZE -DWORDSIZE=32
NAUTYOBJ = $(NAUTY)/gtoolsW.o $(NAUTY)/nautyW1.o $(NAUTY)/nautilW1.o \
           $(NAUTY)/naugraphW1.o $(NAUTY)/schreierW.o $(NAUTY)/naurng.o
SGOBJ    = $(NAUTY)/nausparseW.o        # sparsegraph conversions (nauty_to_sg, sg_to_nauty)

all : geng_coplanar geng_coplanar_instr geng_direct geng_coplanar_lr lrtest

PLANAR   = planarity_arena.c arena.c planarity.h

geng_coplanar : $(NAUTY)/geng.c coplanar.c $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

geng_coplanar_instr : $(NAUTY)/geng.c coplanar.c $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DINSTRUMENT \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

# same generator, but with the left-right planarity test (lrplanar.c) instead
# of nauty's Boyer-Myrvold-style tester
geng_coplanar_lr : $(NAUTY)/geng.c coplanar.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DUSE_LR \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

# compares lrplanar against nauty's tester on graphs read from stdin
lrtest : lrtest.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -o $@ lrtest.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

geng_direct : $(NAUTY)/geng.c direct.c $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) \
	  -DPREPRUNE=direct_preprune -DSUMMARY=direct_summary \
	  -o $@ $(NAUTY)/geng.c direct.c planarity_arena.c arena.c $(NAUTYOBJ)

# experimental: per-parent neighbourhood cache in front of the planarity test
# (-DNBCACHE), optionally with greedy minimisation of non-planar certificates
geng_coplanar_lrc : $(NAUTY)/geng.c coplanar.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DUSE_LR -DNBCACHE \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

geng_coplanar_lrcm : $(NAUTY)/geng.c coplanar.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DUSE_LR -DNBCACHE -DNBCACHE_MIN \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

# --- sparsegraph version of the left-right test (no size limit) and its tests ---
# lrtest_sg links nauty's ORIGINAL planarity.c (malloc-based) as the reference,
# not the arena copy, so that graphs of any size can be tested.
lrplanar_sg.o : lrplanar_sg.c lrplanar_sg.h
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -c lrplanar_sg.c

lrtest_sg : lrtest_sg.c lrplanar_sg.o lrplanar.c lrplanar.h $(NAUTY)/planarity.c $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DNAUTY_PLANARITY_H='"$(NAUTY)/planarity.h"' \
	  -o $@ lrtest_sg.c lrplanar_sg.o lrplanar.c $(NAUTY)/planarity.c $(NAUTYOBJ) $(SGOBJ)

biggraphs : biggraphs.c $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -o $@ biggraphs.c $(NAUTYOBJ) $(SGOBJ)

# the generator with the sparsegraph test as backend (regression test of lrplanar_sg
# on the real workload; slower than geng_coplanar_lr because of the per-call conversion)
geng_coplanar_lrsg : $(NAUTY)/geng.c coplanar.c lrplanar_sg.o $(PLANAR) $(NAUTYOBJ)
	$(CC) $(CFLAGS) $(W1) -I. -I$(NAUTY) -DUSE_LRSG \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar_sg.o planarity_arena.c arena.c $(NAUTYOBJ) $(SGOBJ)

# sanitizer build of the sparsegraph test harness
lrtest_sg_asan : lrtest_sg.c lrplanar_sg.c lrplanar_sg.h lrplanar.c lrplanar.h $(NAUTY)/planarity.c $(NAUTYOBJ)
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(W1) -I. -I$(NAUTY) \
	  -DNAUTY_PLANARITY_H='"$(NAUTY)/planarity.h"' \
	  -o $@ lrtest_sg.c lrplanar_sg.c lrplanar.c $(NAUTY)/planarity.c $(NAUTYOBJ) $(SGOBJ)

tests : lrtest_sg biggraphs geng_coplanar_lrsg lrtest geng_coplanar_lr
	./run_tests.sh

# sanitizer build for testing the hooks and the planarity test (not in 'all')
geng_coplanar_lr_asan : $(NAUTY)/geng.c coplanar.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(W1) -I. -I$(NAUTY) -DUSE_LR \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

clean :
	rm -f geng_coplanar geng_coplanar_instr geng_direct geng_coplanar_lr lrtest geng_coplanar_lr_asan geng_coplanar_lrc geng_coplanar_lrcm \
	  lrplanar_sg.o lrtest_sg biggraphs geng_coplanar_lrsg lrtest_sg_asan

.PHONY : all clean tests
