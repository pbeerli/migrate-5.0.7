bin_PROGRAMS = migrate-n migrate-n-mpi

# Set CFLAGS conditionally for debug
if DEBUG_MODE
AM_CFLAGS = -g -Wall -DDEBUG
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
    mittag_leffler_interpol_data.c romberg.c tables.c

# Shared flags
COMMON_CPPFLAGS = \
    -I../lib/haru/include -I../lib/SFMT-src-1.5.1 \
    -DPRETTY -DNEWVERSION -D_REENTRANT \
    $(ZLIB_CFLAGS) $(LIBHPDF_CFLAGS) \
    -DMERSENNE_TWISTER -DSFMT_MEXP=19937

COMMON_LDADD = \
    ../lib/haru/src/libhpdf.a \
    ../lib/SFMT-src-1.5.1/SFMT.o \
    $(ZLIB_LIBS) $(LIBHPDF_LIBS)

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
    ../lib/SFMT-src-1.5.1/SFMT.o \
    ../lib/haru/src/libhpdf.a
