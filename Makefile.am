SUBDIRS = src

# These are files that are not automatically built, but should be included in dist tarball.
# The SFMT entries were previously just SFMT.c/SFMT.h, but SFMT.h itself
# #includes SFMT-common.h and one of SFMT-params.h/SFMT-params<MEXP>.h/
# SFMT-sse2.h/SFMT-neon.h/SFMT-alti.h depending on platform -- listing
# only SFMT.c/SFMT.h left every one of those out of `make dist`'s
# tarball, breaking the build from a fresh extraction (caught 2026-09-11
# by actually building from one before shipping it). Globbing the whole
# directory's headers (matching the pattern the haru entries below
# already use) covers every platform variant, not just the one this
# machine happens to need.
EXTRA_DIST = \
  lib/SFMT-src-1.5.1/*.c \
  lib/SFMT-src-1.5.1/*.h \
  lib/SFMT-src-1.5.1/LICENSE.txt \
  lib/SFMT-src-1.5.1/README.txt \
  lib/haru/src/*.c \
  lib/haru/src/*.h \
  lib/haru/include/*.h \
  lib/haru/include/hpdf_config.h.in \
  LICENSE \
  HISTORY \
  GUIDANCE \
  README.md \
  example/*

DISTCLEANFILES = \
  lib/SFMT-src-1.5.1/SFMT.o

if !USE_SYSTEM_LIBHPDF
DISTCLEANFILES += \
  lib/haru/src/build/*.o \
  lib/haru/src/libhpdf.a
endif



# Tell Automake we have additional build artifacts
BUILT_SOURCES = \
  lib/SFMT-src-1.5.1/SFMT.o

if !USE_SYSTEM_LIBHPDF
BUILT_SOURCES += lib/haru/src/libhpdf.a
endif

# Hook: build libhpdf.a
#
# The prerequisite list matters: without it, this rule has no way to know
# libhpdf.a is stale once it already exists on disk, so `make` treats any
# already-built libhpdf.a as permanently up to date -- even after editing
# lib/haru/src/*.c -- and silently keeps linking the old one. Found in
# migrate-codex-7 (2026-09-14, ported here): the same missing-dependency
# bug there meant every Autotools-built migrate-n was silently missing
# real libharu fixes made since this .a was first built (June 2025); this
# tree's own libhpdf.a carries the exact same stale June 2025 build.
if !USE_SYSTEM_LIBHPDF
lib/haru/src/libhpdf.a: $(wildcard lib/haru/src/*.c lib/haru/include/*.h)
	@echo "Compiling static libharu (no JPEG, no PNG)..."
	mkdir -p lib/haru/src/build
	cd lib/haru/src && \
	for f in *.c; do \
	  case $$f in \
	    *png*|*jpeg*) continue ;; \
	    *) $(CC) -O3 -I../include -DNOJPEG -DNOPNG -DPRETTY -fPIC -c $$f -o build/$$f.o ;; \
	  esac; \
	done && \
	ar rcs libhpdf.a build/*.o
endif
# Hook: build SFMT.o
lib/SFMT-src-1.5.1/SFMT.o:
	@echo "Compiling SFMT..."
	$(CC) -DMERSENNE_TWISTER -DSFMT_MEXP=19937 -D_REENTRANT -O3 -fPIC -c lib/SFMT-src-1.5.1/SFMT.c -o lib/SFMT-src-1.5.1/SFMT.o

clean-local:
	rm -rf lib/haru/src/build
