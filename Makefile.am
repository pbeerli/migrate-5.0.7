SUBDIRS = src

# These are files that are not automatically built, but should be included in dist tarball
EXTRA_DIST = \
  lib/SFMT-src-1.5.1/SFMT.c \
  lib/SFMT-src-1.5.1/SFMT.h \
  lib/haru/src/*.c \
  lib/haru/include/*.h

DISTCLEANFILES = \
  lib/SFMT-src-1.5.1/SFMT.o \
  lib/haru/src/build/*.o \
  lib/haru/src/libhpdf.a

# Tell Automake we have additional build artifacts
BUILT_SOURCES = \
  lib/SFMT-src-1.5.1/SFMT.o \
  lib/haru/src/libhpdf.a

# Hook: build libhpdf.a
lib/haru/src/libhpdf.a:
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

# Hook: build SFMT.o
lib/SFMT-src-1.5.1/SFMT.o:
	@echo "Compiling SFMT..."
	$(CC) -DMERSENNE_TWISTER -DSFMT_MEXP=19937 -D_REENTRANT -O3 -fPIC -c lib/SFMT-src-1.5.1/SFMT.c -o lib/SFMT-src-1.5.1/SFMT.o

clean-local:
	rm -rf lib/haru/src/build
