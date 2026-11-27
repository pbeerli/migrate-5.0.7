/* Checks that a libharu writes a valid /ToUnicode CMap for a Unicode (UTF-8,
   Identity-H) TrueType font, the fix that libharu-pdfua carries in
   vendor/libharu/src/hpdf_font_cid.c (see its NOTICE.md). Stock libharu
   2.4.5 writes the /ToUnicode stream with cidrange/cidchar operators, which
   are not valid there (PDF32000-1:2008 9.10.3 requires bfrange/bfchar):
   PDF readers and veraPDF then cannot map migrate's symbol glyphs (Theta,
   Greek letters) to Unicode, and the report is not PDF/UA-1.

   Used by configure (--with-system-libhpdf) and CMake
   (MIGRATE_USE_SYSTEM_HARU) before they link a system libharu. Build:
     cc -I<hpdf include> -Isrc tools/check_libhpdf_tounicode.c -lhpdf -lz -lm
   Exit status: 0 the fix is present, 1 it is missing, 2 the test failed to
   run. The other libharu-pdfua fixes (encoder buffer bound, PDF/A enum,
   deflate buffer size) do not change a correct document and cannot be
   tested this way. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hpdf.h"
#include "pdf_fonts_embedded.c"   /* migrate's embedded Liberation Sans */

static void error_handler (HPDF_STATUS error_no, HPDF_STATUS detail_no, void *user_data)
{
  (void) user_data;
  fprintf (stderr, "libharu error %04X (detail %u)\n", (unsigned) error_no, (unsigned) detail_no);
  exit (2);
}

int main (void)
{
  HPDF_Doc doc = HPDF_New (error_handler, NULL);
  HPDF_Page page;
  HPDF_Font font;
  const char *name;
  HPDF_UINT32 size;
  HPDF_BYTE *buf;
  const char *bf, *cid;
  if (doc == NULL)
    return 2;
  HPDF_SetCompressionMode (doc, HPDF_COMP_NONE);
  name = HPDF_LoadTTFontFromMemory (doc, pdf_font_helv_ttf, pdf_font_helv_ttf_len, HPDF_TRUE);
  if (name == NULL || HPDF_UseUTFEncodings (doc) != HPDF_OK)
    return 2;
  font = HPDF_GetFont (doc, name, "UTF-8");
  page = HPDF_AddPage (doc);
  HPDF_Page_SetFontAndSize (page, font, 12);
  HPDF_Page_BeginText (page);
  HPDF_Page_TextOut (page, 50, 700, "\xce\x98 = 0.01");   /* U+0398 Theta */
  HPDF_Page_EndText (page);
  if (HPDF_SaveToStream (doc) != HPDF_OK)
    return 2;
  size = HPDF_GetStreamSize (doc);
  buf = (HPDF_BYTE *) malloc (size + 1);
  if (buf == NULL)
    return 2;
  {
    HPDF_STATUS st = HPDF_ReadFromStream (doc, buf, &size);
    if (st != HPDF_OK && st != HPDF_STREAM_EOF)
      return 2;
  }
  buf[size] = 0;
  /* the PDF holds NUL bytes (embedded font): search the whole buffer */
  bf = cid = NULL;
  {
    HPDF_UINT32 i;
    for (i = 0; i + 13 <= size && (bf == NULL || cid == NULL); i++)
      {
        if (bf == NULL && (memcmp (buf + i, "beginbfrange", 12) == 0 || memcmp (buf + i, "beginbfchar", 11) == 0))
          bf = (const char *) buf + i;
        if (cid == NULL && (memcmp (buf + i, "begincidrange", 13) == 0 || memcmp (buf + i, "begincidchar", 12) == 0))
          cid = (const char *) buf + i;
      }
  }
  HPDF_Free (doc);
  if (bf != NULL && cid == NULL)
    {
      printf ("libharu writes a valid /ToUnicode CMap (bfrange/bfchar)\n");
      return 0;
    }
  printf ("libharu lacks the /ToUnicode fix: %s\n",
          cid != NULL ? "the CMap uses cidrange/cidchar" : "no bfrange/bfchar CMap found");
  return 1;
}
