AllExecs: AllLibs

CFLAGS ?= -O3 -march=native -mieee-fp -pipe -Wall -Wcast-align

ifdef CT_OUTPUT
  CFLAGS += -DCT_OUTPUT
endif

ifdef OBSREFF_OUTPUT
  CFLAGS += -DOBSREFF_OUTPUT
endif

ifdef DUAL_PINF
  CFLAGS += -DDUAL_PINF
endif

include src/lib/Makefile.in src/bin/Makefile.in src/doc/Makefile.in
