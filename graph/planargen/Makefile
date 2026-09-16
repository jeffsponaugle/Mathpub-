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

# sanitizer build for testing the hooks and the planarity test (not in 'all')
geng_coplanar_lr_asan : $(NAUTY)/geng.c coplanar.c lrplanar.c lrplanar.h $(PLANAR) $(NAUTYOBJ)
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer $(W1) -I. -I$(NAUTY) -DUSE_LR \
	  -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune -DSUMMARY=coplanar_summary \
	  -o $@ $(NAUTY)/geng.c coplanar.c lrplanar.c planarity_arena.c arena.c $(NAUTYOBJ)

clean :
	rm -f geng_coplanar geng_coplanar_instr geng_direct geng_coplanar_lr lrtest geng_coplanar_lr_asan

.PHONY : all clean
