#!/bin/sh
# Syncs the PDF/UA-1 tagging module (hpdf_ua_context.c/hpdf_ua_metadata.c/
# hpdf_ua_docmeta.c/hpdf_ua_structure.c/hpdf_ua_stack.c/hpdf_ua.h/hpdf_ua_private.h) from the
# public libharu-pdfua project (https://github.com/pbeerli/libharu-pdfua)
# into src/, and (optionally) libharu's own hpdf_font_cid.c ToUnicode fix in
# the other direction, into lib/haru/.
#
# Why this exists: this module isn't part of libharu itself (see
# hpdf_ua_context.c's own file comment -- it only calls libharu's public
# API, no vendored libharu header is touched to make it work), so it isn't
# vendored under lib/haru/ alongside plain libharu; historically it was
# copy-pasted flat into src/ once, by hand, with no record of which
# libharu-pdfua commit it came from and no repeatable way to re-pull a
# later improvement (e.g. HPDF_UA_TagAnnotation(), added to libharu-pdfua
# after this project's copy was made). This script makes that a real,
# repeatable, provenance-tracked step instead, mirroring the existing
# tools/generate_embedded_fonts.sh precedent (checked-in output, no
# network/build-time dependency -- migrate targets HPC clusters that may
# build offline, so this deliberately stays a vendor-and-check-in model,
# not a git submodule or CMake FetchContent).
#
# The only difference between migrate's copies and libharu-pdfua's own
# files is two include lines (migrate keeps this module flat in src/,
# libharu-pdfua uses its own include/hpdf_ua/ + src/ua/ layout) -- this
# script applies that rewrite mechanically, so migrate's copies should
# never need hand-editing: change libharu-pdfua, re-run this script.
#
# Usage:
#   tools/sync_libharu_pdfua.sh [path-to-libharu-pdfua-checkout]
# Defaults to $LIBHARU_PDFUA_DIR or ~/src/libharu-pdfua.
set -e
cd "$(dirname "$0")/.."

SRC="${1:-${LIBHARU_PDFUA_DIR:-$HOME/src/libharu-pdfua}}"
if [ ! -d "$SRC/.git" ]; then
    echo "error: $SRC is not a libharu-pdfua git checkout" >&2
    exit 1
fi

SRC_COMMIT=$(cd "$SRC" && git rev-parse --short=12 HEAD)
SRC_DIRTY=""
if ! (cd "$SRC" && git diff --quiet && git diff --cached --quiet); then
    SRC_DIRTY=" (with local uncommitted changes -- not a clean sync point)"
fi

NOTE="/* Synced from libharu-pdfua commit $SRC_COMMIT via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */"

# $1 = source path relative to $SRC, $2 = dest path relative to repo root,
# $3 = sed expression to rewrite this project's one include-path difference
# (empty = copy verbatim).
sync_file () {
    src_file="$SRC/$1"
    dest_file="$2"
    rewrite="$3"

    if [ ! -f "$src_file" ]; then
        echo "error: missing $src_file" >&2
        exit 1
    fi

    {
        echo "$NOTE"
        if [ -n "$rewrite" ]; then
            sed "$rewrite" "$src_file"
        else
            cat "$src_file"
        fi
    } > "$dest_file"
}

sync_file src/ua/hpdf_ua_context.c   src/hpdf_ua_context.c
sync_file src/ua/hpdf_ua_metadata.c  src/hpdf_ua_metadata.c
sync_file src/ua/hpdf_ua_structure.c src/hpdf_ua_structure.c
sync_file src/ua/hpdf_ua_stack.c     src/hpdf_ua_stack.c
sync_file src/ua/hpdf_ua_docmeta.c   src/hpdf_ua_docmeta.c \
    's|#include "hpdf_ua/hpdf_ua.h"|#include "hpdf_ua.h"|'
sync_file src/ua/hpdf_ua_private.h   src/hpdf_ua_private.h \
    's|#include "hpdf_ua/hpdf_ua.h"|#include "hpdf_ua.h"|'
sync_file include/hpdf_ua/hpdf_ua.h  src/hpdf_ua.h

echo "Synced from libharu-pdfua commit $SRC_COMMIT$SRC_DIRTY."
echo "Review with: git diff -- src/hpdf_ua_context.c src/hpdf_ua_metadata.c src/hpdf_ua_docmeta.c src/hpdf_ua_structure.c src/hpdf_ua_stack.c src/hpdf_ua.h src/hpdf_ua_private.h"
echo "Then rebuild and re-run the full verification pass (CTest, sanitizer build, veraPDF on real examples) before committing."
