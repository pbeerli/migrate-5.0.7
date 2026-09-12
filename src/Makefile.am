bin_PROGRAMS = migrate-n migrate-n-mpi

# Set CFLAGS conditionally for debug
if DEBUG_MODE
AM_CFLAGS = -g -Wall -DDEBUG -fsanitize=address -fno-omit-frame-pointer
else
AM_CFLAGS = -O3 -DNDEBUG
endif

COMMON_SOURCES = \
    assignment.c autotune.c bayes.c calculator.c correlation.c data.c \
    growth.c haplotype.c hash.c heating.c histogram.c kernel_smooth.c laguerre.c main.c marginallike.c \
    mcmc1.c mcmc2.c menu.c migevents.c migrate_mpi.c mutationmodel.c \
    options.c pretty.c priors.c random.c reporter.c sequence.c seqerror.c \
    savitzky_golay.c \
    sighandler.c skyline.c skyparam.c slice.c sort.c speciate.c tools.c \
    tree.c uep.c watterson.c world.c znzlib.c mittag_leffler.c hermite_interpoly.c \
    mittag_leffler_interpol_data.c romberg.c tables.c mlalpha.c windowupdate.c

# Not compiled directly, but every one of them is #included by the
# sources above -- listing them here is what makes `make dist` actually
# package them. Without this, "make dist"'s tarball builds nothing at
# all (every .c's first #include fails) despite `make dist` itself
# reporting success; caught 2026-09-11 by extracting a fresh tarball
# into a clean directory and building from it before shipping it.
noinst_HEADERS = \
    assignment.h autotune.h bayes.h calculator.h conf.h correlation.h \
    data.h definitions.h growth.h haplotype.h hash.h heating.h \
    hermite_interpoly.h histogram.h kernel_smooth.h laguerre.h \
    marginallike.h mcmc.h mcmc2.h menu.h migevents.h migrate_mpi.h \
    migration.h mittag_leffler_interpol_data.h mittag_leffler.h mlalpha.h \
    mutationmodel.h options.h pretty.h priors.h random.h reporter.h \
    romberg.h savitzky_golay.h seqerror.h sequence.h sighandler.h \
    skyline.h skyparam.h slice.h sort.h speciate.h tables.h tools.h \
    tree.h uep.h watterson.h windows_timveval.h windowupdate.h world.h \
    znzlib.h

# Shared flags (common defines + local includes)
COMMON_CPPFLAGS = \
    -I../lib/SFMT-src-1.5.1 \
    -DPRETTY -DNEWVERSION -D_REENTRANT \
    -DMERSENNE_TWISTER -DSFMT_MEXP=19937

COMMON_LDADD = \
    ../lib/SFMT-src-1.5.1/SFMT.o

# --- zlib selection (optional, you already have USE_SYSTEM_ZLIB) ---
if USE_SYSTEM_ZLIB
COMMON_CPPFLAGS += $(ZLIB_CFLAGS)
else
# If you vendor zlib, add your bundled zlib include/lib here
# COMMON_CPPFLAGS += -I../lib/zlib/include
# COMMON_LDADD    += ../lib/zlib/libz.a
endif

# --- libharu selection ---
if USE_SYSTEM_LIBHPDF
COMMON_CPPFLAGS += $(LIBHPDF_CFLAGS)
COMMON_LDADD    += $(LIBHPDF_LIBS)
else
COMMON_CPPFLAGS += -I../lib/haru/include -I../lib/haru/build/include
COMMON_LDADD    += ../lib/haru/src/libhpdf.a
endif
COMMON_LDADD    += $(ZLIB_LIBS)

# migrate-n
migrate_n_SOURCES = $(COMMON_SOURCES)
migrate_n_CPPFLAGS = $(COMMON_CPPFLAGS)
migrate_n_CFLAGS = $(AM_CFLAGS)
migrate_n_LDADD = $(COMMON_LDADD)

# migrate-n-mpi (parallel)
migrate_n_mpi_SOURCES = $(COMMON_SOURCES)
migrate_n_mpi_CPPFLAGS = $(COMMON_CPPFLAGS) $(MPICFLAGS) \
  -DMPI -DSLOWNET -DNEWREPLICANT -DMPIREPLICANT -DMPIDATAONDEMAND -DIPROBE
migrate_n_mpi_CFLAGS = $(AM_CFLAGS)
migrate_n_mpi_LDADD = $(COMMON_LDADD) $(MPILIBS)
migrate_n_mpi_CC = $(MPICC)

migrate_n_LDADD += $(MATHLIBS)
migrate_n_mpi_LDADD += $(MATHLIBS)

# Trigger build of bundled sources
BUILT_SOURCES = \
    ../lib/SFMT-src-1.5.1/SFMT.o

if !USE_SYSTEM_LIBHPDF
BUILT_SOURCES += ../lib/haru/src/libhpdf.a
endif
