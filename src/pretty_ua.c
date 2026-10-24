#include "pretty_ua.h"

/* See pretty_ua.h for the design rationale. A thin wrapper around
 * libharu-pdfua's HPDF_UA_Stack (hpdf_ua_stack.c), which holds the open
 * elements, the open content/artifact span and the two independent
 * suspend/resume states (page break, nested element on the same page).
 * Module-global, single-threaded state, matching pretty.c's own existing
 * globals (`doc`, `page`, `helv`, ...) -- the PDF report is only ever
 * generated once, on the master rank, never concurrently.
 *
 * What this wrapper adds is tolerance: pretty.c's error_handler()
 * longjmps out of the whole PDF output on any libharu error, and the
 * stack raises one for misuse (a second begin while a span is open, an
 * end or pop with nothing open). The old glue silently ignored such
 * calls, and its callers rely on that (e.g. pretty_ua_end_content()
 * after a pop, which now already closed the span), so every call checks
 * the stack state first and quietly does nothing when it would be
 * invalid. */

static HPDF_UA_Context g_ctx = NULL;
static HPDF_UA_Stack g_stack = NULL;

void
pretty_ua_start_document (HPDF_Doc pdf, const char *title)
{
  HPDF_SetInfoAttr (pdf, HPDF_INFO_TITLE, title);
  HPDF_UA_SetDocumentLanguage (pdf, "en-US");
  HPDF_UA_SetDisplayDocTitle (pdf, HPDF_TRUE);
  HPDF_UA_AddMetadata (pdf);

  g_ctx = HPDF_UA_NewContext (pdf);
  if (!g_ctx)
    return; /* degrade to "no tagging" rather than crash -- pretty.c's
              * drawing calls are unaffected either way. */
  g_stack = HPDF_UA_NewStack (g_ctx);
  if (!g_stack)
    {
      HPDF_UA_FreeContext (g_ctx);
      g_ctx = NULL;
      return;
    }
  HPDF_UA_Stack_Push (g_stack, HPDF_UA_ROLE_DOCUMENT);
}

void
pretty_ua_finish_document (void)
{
  if (g_stack)
    {
      HPDF_UA_FreeStack (g_stack); /* closes open spans, ends all elements */
      g_stack = NULL;
    }
  if (g_ctx)
    {
      HPDF_UA_FreeContext (g_ctx);
      g_ctx = NULL;
    }
}

HPDF_UA_StructElem
pretty_ua_push (HPDF_UA_StructType role)
{
  if (!g_stack)
    return NULL;
  return HPDF_UA_Stack_Push (g_stack, role);
}

void
pretty_ua_pop (void)
{
  if (!g_stack || HPDF_UA_Stack_Depth (g_stack) == 0)
    return;
  HPDF_UA_Stack_Pop (g_stack); /* also ends an open content span */
}

void
pretty_ua_set_alt_text (const char *alt_text)
{
  if (!g_stack || HPDF_UA_Stack_Depth (g_stack) == 0)
    return;
  HPDF_UA_Stack_SetAlternateText (g_stack, alt_text);
}

/// Same idiom as pretty_ua_set_alt_text() above (applies to whichever
/// element is currently on top of the stack), but sets /ActualText
/// instead of /Alt -- for a TD/TH whose visible content was drawn
/// through the Greek/math "symb" (Identity-H/CID) font, giving checkers
/// a direct Unicode text equivalent that doesn't depend on how well they
/// walk that font's own content-stream text-showing operators. Pass the
/// same UTF-8 text the glyph itself represents (e.g. "Theta_1", the
/// caller's own real label -- not the raw glyph string), see call sites
/// in pretty.c/pretty_table.c.
void
pretty_ua_set_actual_text (const char *actual_text)
{
  if (!g_stack || HPDF_UA_Stack_Depth (g_stack) == 0)
    return;
  HPDF_UA_Stack_SetActualText (g_stack, actual_text);
}

void
pretty_ua_set_table_header_scope (HPDF_UA_TableScope scope)
{
  if (!g_stack || HPDF_UA_Stack_Depth (g_stack) == 0)
    return;
  HPDF_UA_Stack_SetTableHeaderScope (g_stack, scope);
}

/* a content span needs an element to belong to and no other open span */
static HPDF_BOOL
pretty_ua_can_begin_span (void)
{
  return g_stack && HPDF_UA_Stack_Depth (g_stack) > 0
         && !HPDF_UA_Stack_ContentOpen (g_stack)
         && !HPDF_UA_Stack_ArtifactOpen (g_stack);
}

void
pretty_ua_begin_content (HPDF_Page page)
{
  if (pretty_ua_can_begin_span ())
    HPDF_UA_Stack_BeginContent (g_stack, page);
}

void
pretty_ua_begin_content_with_actual_text (HPDF_Page page, const char *actual_text)
{
  /* a NULL actual_text behaves like plain pretty_ua_begin_content() */
  if (pretty_ua_can_begin_span ())
    HPDF_UA_Stack_BeginContentWithActualText (g_stack, page, actual_text);
}

void
pretty_ua_end_content (HPDF_Page page)
{
  (void) page; /* the stack remembers the span's page */
  if (g_stack && HPDF_UA_Stack_ContentOpen (g_stack))
    HPDF_UA_Stack_EndContent (g_stack);
}

void
pretty_ua_begin_artifact (HPDF_Page page)
{
  if (g_stack && !HPDF_UA_Stack_ArtifactOpen (g_stack)
      && !HPDF_UA_Stack_ContentOpen (g_stack))
    HPDF_UA_Stack_BeginArtifact (g_stack, page, HPDF_UA_ARTIFACT_NONE);
}

void
pretty_ua_end_artifact (HPDF_Page page)
{
  (void) page;
  if (g_stack && HPDF_UA_Stack_ArtifactOpen (g_stack))
    HPDF_UA_Stack_EndArtifact (g_stack);
}

/* The suspend/resume pairs map one to one onto the stack's: the content
 * element is popped for good and a new sibling of the same role is pushed
 * on resume (an element's content lives on one page only), and the two
 * pairs keep separate saved state because their windows nest (see
 * hpdf_ua.h, HPDF_UA_Stack_SuspendContent()). */
void
pretty_ua_suspend_for_new_page (HPDF_Page old_page)
{
  (void) old_page;
  if (g_stack)
    HPDF_UA_Stack_SuspendForNewPage (g_stack);
}

void
pretty_ua_resume_after_new_page (HPDF_Page new_page)
{
  if (g_stack)
    HPDF_UA_Stack_ResumeAfterNewPage (g_stack, new_page);
}

void
pretty_ua_suspend_content_same_page (HPDF_Page page)
{
  (void) page;
  if (g_stack)
    HPDF_UA_Stack_SuspendContent (g_stack);
}

void
pretty_ua_resume_content_same_page (HPDF_Page page)
{
  /* `page` is whatever pretty.c's global `page` holds now -- possibly a
   * later page than at the suspend, if the nested drawing in between
   * broke the page; the resumed sibling goes on the current page. */
  if (g_stack)
    HPDF_UA_Stack_ResumeContent (g_stack, page);
}
