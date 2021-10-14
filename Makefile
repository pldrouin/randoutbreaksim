AllExecs: AllLibs

CFLAGS ?= -O3 -march=native -mieee-fp -pipe -Wall -Wcast-align

ifdef NUMEVENTSSTATS
  CFLAGS += -DNUMEVENTSSTATS
endif

ifdef CT_OUTPUT
  CFLAGS += -DCT_OUTPUT
endif

ifdef OBSREFF_OUTPUT
  CFLAGS += -DOBSREFF_OUTPUT
endif

ifdef DUAL_PINF
  CFLAGS += -DDUAL_PINF

  ifdef SEC_INF_TIMELINES
    CFLAGS += -DSEC_INF_TIMELINES
  endif
endif

include src/lib/Makefile.in src/bin/Makefile.in src/doc/Makefile.in
