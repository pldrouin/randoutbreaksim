AllExecs: AllLibs

CFLAGS ?= -O3 -march=native -mieee-fp -pipe -Wall -Wcast-align

ifdef CT_OUTPUT
  CFLAGS += -DCT_OUTPUT
endif

include src/lib/Makefile.in src/bin/Makefile.in src/doc/Makefile.in
