#include "pretty_ua.h"

/* See pretty_ua.h for the design rationale. Module-global, single-
 * threaded state, matching pretty.c's own existing globals (`doc`, `page`,
 * `helv`, ...) -- the PDF report is only ever generated once, on the
 * master rank, never concurrently. */

#define PRETTY_UA_MAX_DEPTH 64

static HPDF_UA_Context g_ctx = NULL;
static HPDF_UA_StructElem g_stack[PRETTY_UA_MAX_DEPTH];
static HPDF_UA_StructType g_stack_role[PRETTY_UA_MAX_DEPTH];
static int g_depth = 0;
static int g_content_open = 0;  /* guards against a stray unmatched end_content() */
static int g_artifact_open = 0;

/* Two INDEPENDENT sets of "suspended" state, one per public suspend/
 * resume pair (see pretty_ua.h). These must NOT share storage: a
 * pretty_ua_suspend_content_same_page()/pretty_ua_resume_content_same_page()
 * pair, used around one nested table's whole rendering, routinely spans
 * zero or more pdf_new_page() calls (and hence zero or more
 * pretty_ua_suspend_for_new_page()/pretty_ua_resume_after_new_page()
 * pairs) triggered by that table's own pdf_advance() calls -- so the two
 * pairs' "currently outstanding suspend" windows genuinely nest, they do
 * not merely alternate. An earlier version of this file used one shared
 * set of globals for both, reasoning that the inner (page-break) pair
 * would always find nothing to suspend of its own (since a nested
 * table's TD/TH content is always closed before it ever calls
 * pdf_advance()) and so would leave the shared state alone. That
 * reasoning missed that pretty_ua_suspend_for_new_page() *unconditionally*
 * overwrites the shared "what to resume" slot with its own (empty)
 * findings every time pdf_new_page() runs, even when there is nothing
 * for it to actually suspend -- silently clobbering the outer pair's
 * still-pending state whenever a page break happened to occur somewhere
 * in between. The result: whatever the outer pair was supposed to
 * reopen on resume silently never was, leaving real content genuinely
 * untagged -- not a crash this time, but a real regression, caught by
 * this project's own veraPDF-delta verification (failed checks went up,
 * not down, on a case where the affected table happens to straddle a
 * page break) before being committed. Two separate sets of globals make
 * the two windows independent, however they nest. */
static int g_suspended_content = 0;
static HPDF_UA_StructType g_suspended_content_role;
static int g_suspended_artifact = 0;

static int g_nested_suspended_content = 0;
static HPDF_UA_StructType g_nested_suspended_content_role;
static int g_nested_suspended_artifact = 0;

void
pretty_ua_start_document (HPDF_Doc pdf, const char *title)
{
  HPDF_UA_StructElem root;

  HPDF_SetInfoAttr (pdf, HPDF_INFO_TITLE, title);
  HPDF_UA_SetDocumentLanguage (pdf, "en-US");
  HPDF_UA_SetDisplayDocTitle (pdf, HPDF_TRUE);
  HPDF_UA_AddMetadata (pdf);

  g_ctx = HPDF_UA_NewContext (pdf);
  g_depth = 0;
  g_content_open = 0;
  g_artifact_open = 0;

  if (!g_ctx)
    return; /* degrade to "no tagging" rather than crash -- pretty.c's
              * drawing calls are unaffected either way. */

  root = HPDF_UA_BeginStructureElement (g_ctx, NULL, HPDF_UA_ROLE_DOCUMENT);
  if (root)
    {
      g_stack_role[g_depth] = HPDF_UA_ROLE_DOCUMENT;
      g_stack[g_depth++] = root;
    }
}

void
pretty_ua_finish_document (void)
{
  while (g_depth > 0)
    {
      g_depth--;
      HPDF_UA_EndStructureElement (g_ctx, g_stack[g_depth]);
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
  HPDF_UA_StructElem parent;
  HPDF_UA_StructElem elem;

  if (!g_ctx)
    return NULL;

  parent = (g_depth > 0) ? g_stack[g_depth - 1] : NULL;
  elem = HPDF_UA_BeginStructureElement (g_ctx, parent, role);
  if (elem && g_depth < PRETTY_UA_MAX_DEPTH)
    {
      g_stack_role[g_depth] = role;
      g_stack[g_depth++] = elem;
    }

  return elem;
}

void
pretty_ua_pop (void)
{
  if (!g_ctx || g_depth == 0)
    return;

  g_depth--;
  HPDF_UA_EndStructureElement (g_ctx, g_stack[g_depth]);
}

void
pretty_ua_set_alt_text (const char *alt_text)
{
  if (!g_ctx || g_depth == 0)
    return;
  HPDF_UA_SetAlternateText (g_ctx, g_stack[g_depth - 1], alt_text);
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
  if (!g_ctx || g_depth == 0)
    return;
  HPDF_UA_SetActualText (g_ctx, g_stack[g_depth - 1], actual_text);
}

void
pretty_ua_set_table_header_scope (HPDF_UA_TableScope scope)
{
  if (!g_ctx || g_depth == 0)
    return;
  HPDF_UA_SetTableHeaderScope (g_ctx, g_stack[g_depth - 1], scope);
}

void
pretty_ua_begin_content (HPDF_Page page)
{
  if (!g_ctx || g_depth == 0 || g_content_open)
    return;
  if (HPDF_UA_BeginMarkedContent (g_ctx, page, g_stack[g_depth - 1]) == HPDF_OK)
    g_content_open = 1;
}

void
pretty_ua_begin_content_with_actual_text (HPDF_Page page, const char *actual_text)
{
  if (!actual_text)
    {
      /* No text given -- behave exactly like plain pretty_ua_begin_content()
         rather than failing (HPDF_UA_BeginMarkedContentWithActualText()
         itself rejects a NULL actual_text). */
      pretty_ua_begin_content (page);
      return;
    }
  if (!g_ctx || g_depth == 0 || g_content_open)
    return;
  if (HPDF_UA_BeginMarkedContentWithActualText (g_ctx, page, g_stack[g_depth - 1],
                                                 actual_text) == HPDF_OK)
    g_content_open = 1;
}

void
pretty_ua_end_content (HPDF_Page page)
{
  if (!g_content_open)
    return;
  HPDF_UA_EndMarkedContent (g_ctx, page);
  g_content_open = 0;
}

void
pretty_ua_begin_artifact (HPDF_Page page)
{
  if (!g_ctx || g_artifact_open)
    return;
  if (HPDF_UA_BeginArtifact (g_ctx, page) == HPDF_OK)
    g_artifact_open = 1;
}

void
pretty_ua_end_artifact (HPDF_Page page)
{
  if (!g_artifact_open)
    return;
  HPDF_UA_EndArtifact (g_ctx, page);
  g_artifact_open = 0;
}

/* Shared engine behind both public suspend/resume pairs (see
 * pretty_ua.h), parameterized on which pair's own (independent -- see
 * the file header comment above) state slots to use. Whatever is open
 * gets closed for good -- the struct element itself popped, not merely
 * suspended, because HPDF_UA_BeginMarkedContent() hard-enforces "one
 * struct element's content lives on only one page" (see
 * hpdf_ua_structure.c's file header) and there is no way to know in
 * advance, at suspend time, whether the gap until resume will turn out
 * to cross a page or not. Always popping and, on resume, always pushing
 * a brand-new sibling of the same role sidesteps the question entirely:
 * it is correct whether or not a page break happened in between, and
 * however many happened. */
static void
pretty_ua_suspend_and_pop (HPDF_Page page, int *suspended_content,
                           HPDF_UA_StructType *suspended_content_role,
                           int *suspended_artifact)
{
  *suspended_content  = g_content_open;
  *suspended_artifact = g_artifact_open;

  if (*suspended_content)
    {
      *suspended_content_role = g_stack_role[g_depth - 1];
      pretty_ua_end_content (page);
      pretty_ua_pop ();
    }

  if (*suspended_artifact)
    pretty_ua_end_artifact (page);
}

static void
pretty_ua_resume_as_new_sibling (HPDF_Page page, int *suspended_content,
                                  HPDF_UA_StructType suspended_content_role,
                                  int *suspended_artifact)
{
  if (*suspended_content)
    {
      pretty_ua_push (suspended_content_role);
      pretty_ua_begin_content (page);
      *suspended_content = 0;
    }

  if (*suspended_artifact)
    {
      pretty_ua_begin_artifact (page);
      *suspended_artifact = 0;
    }
}

void
pretty_ua_suspend_for_new_page (HPDF_Page old_page)
{
  pretty_ua_suspend_and_pop (old_page, &g_suspended_content,
                             &g_suspended_content_role, &g_suspended_artifact);
  /* Both spans are now fully closed -- pdf_new_page() is free to open
   * and close its own border/header/footer Artifact span on the new
   * page with no risk of overlapping a still-open span from before. */
}

void
pretty_ua_resume_after_new_page (HPDF_Page new_page)
{
  pretty_ua_resume_as_new_sibling (new_page, &g_suspended_content,
                                   g_suspended_content_role, &g_suspended_artifact);
}

void
pretty_ua_suspend_content_same_page (HPDF_Page page)
{
  pretty_ua_suspend_and_pop (page, &g_nested_suspended_content,
                             &g_nested_suspended_content_role,
                             &g_nested_suspended_artifact);
}

void
pretty_ua_resume_content_same_page (HPDF_Page page)
{
  /* `page` is simply whatever the caller's own global `page` variable
   * holds *right now* -- possibly a different page than at the matching
   * suspend call, if the nested element(s) drawn in between happened to
   * overflow onto one or more new pages of their own (each already
   * handled correctly and independently by pdf_new_page()'s own
   * suspend/resume pair around its own border artifact -- using its own,
   * separate state slots, see the file header comment above). Either
   * way, resuming here always means a fresh sibling on whatever page is
   * current now -- see pretty_ua_suspend_and_pop()'s comment for why. */
  pretty_ua_resume_as_new_sibling (page, &g_nested_suspended_content,
                                   g_nested_suspended_content_role,
                                   &g_nested_suspended_artifact);
}
