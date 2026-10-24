#ifndef MIGRATE_PRETTY_UA_H
#define MIGRATE_PRETTY_UA_H
/*------------------------------------------------------------------------
  PDF/UA-1 tagging glue for pretty.c (integrating the libharu-pdfua side
  project, see plan.md's "Independent Haru project" section / Milestone 5).

  Design: pretty.c's existing drawing primitives (pdf_printf() and its
  family, pdf_draw_line(), pdf_fill_stroke(), all funneling
  through the single global `page`) are NOT changed to take a structure
  element parameter each -- that would mean touching every one of
  pretty.c's ~90 report-producing functions individually, all at once, a
  large, high-risk change. Instead, this module keeps a small stack (a
  thin wrapper around libharu-pdfua's HPDF_UA_Stack, see pretty_ua.c) of
  "currently open" structure elements; a handful of call sites at each
  function's natural section boundaries (start/end of a table, a row, a
  cell, a figure) push/pop that stack and bracket a block of otherwise
  completely unchanged drawing calls with pretty_ua_begin_content()/
  end_content() (real content) or pretty_ua_begin_artifact()/end_artifact()
  (purely decorative marks -- borders, the migrate logo). Nothing in this
  header changes pretty.c's actual visual output; every currently-unwrapped
  function keeps drawing exactly as it always has -- untagged, same as
  before this integration -- until it too gets a push/pop pair added.

  Safe to call before pretty_ua_start_document() (no-ops) and safe to nest
  pretty_ua_begin_content()/begin_artifact() calls incorrectly (each is
  independently guarded) -- this module is meant to be extended
  incrementally, function by function, without risking already-instrumented
  code.
------------------------------------------------------------------------*/
#include "hpdf_ua.h"

/// Creates the tagging context for `pdf`, sets /Lang, /DisplayDocTitle,
/// and a real XMP /Metadata stream (declaring PDF/UA-1 conformance), sets
/// the document Info /Title to `title`, and pushes an initial Document
/// structure element (so any content tagged before a more specific
/// element is pushed still gets a sensible parent). Call once, right
/// after HPDF_New().
void pretty_ua_start_document (HPDF_Doc pdf, const char *title);

/// Pops any remaining open elements and frees the tagging context. Call
/// once, right before HPDF_SaveToFile().
void pretty_ua_finish_document (void);

/// Pushes a new structure element as a child of whatever is currently on
/// top of the stack (or the Document root if nothing more specific has
/// been pushed), making it the new top. Returns the new element (mostly
/// useful for pretty_ua_set_alt_text()/set_table_header_scope() below);
/// safe to ignore the return value otherwise. Returns NULL (and leaves
/// the stack unchanged) if pretty_ua_start_document() was never called
/// or the underlying allocation failed -- callers do not need to check
/// this themselves, later content simply attaches to whatever the
/// previous top was, or stays untagged if there is none.
HPDF_UA_StructElem pretty_ua_push (HPDF_UA_StructType role);

/// Pops (and ends) the current top of the stack. A no-op if the stack is
/// already empty.
void pretty_ua_pop (void);

/// Sets /Alt on the current top of the stack (typically called right
/// after pretty_ua_push(HPDF_UA_ROLE_FIGURE)). No-op if the stack is
/// empty.
void pretty_ua_set_alt_text (const char *alt_text);

/// Sets /ActualText on the current top of the stack (typically called
/// right after pretty_ua_push(HPDF_UA_ROLE_TD)/(HPDF_UA_ROLE_TH), around
/// content drawn through the Greek/math "symb" CID font) -- a direct
/// Unicode text equivalent for checkers that can't otherwise recover the
/// cell's text. No-op if the stack is empty.
void pretty_ua_set_actual_text (const char *actual_text);

/// Sets /Scope on the current top of the stack (typically called right
/// after pretty_ua_push(HPDF_UA_ROLE_TH)). No-op if the stack is empty.
void pretty_ua_set_table_header_scope (HPDF_UA_TableScope scope);

/// Same as pretty_ua_begin_content() below, but also writes `actual_text`
/// as this span's own /ActualText, directly in the content stream (not
/// just on the structure element the way pretty_ua_set_actual_text()
/// does) -- use this instead of pretty_ua_set_actual_text() wherever a
/// checker needs to recover the text and can't reliably do so by
/// reading the content stream itself, such as content drawn through the
/// "symb" CID font (Theta/mu/alpha): a struct-element-only /ActualText
/// turned out not to be enough in practice for at least one real
/// checker (avalpdf/pdfix-sdk), confirmed by direct testing.
void pretty_ua_begin_content_with_actual_text (HPDF_Page page, const char *actual_text);

/// Opens one marked-content span on `page`, tied to the current top of
/// the stack, wrapping whatever ordinary pretty.c drawing calls
/// (pdf_printf(), pdf_draw_line(), pdf_fill_stroke(), ...) happen between this
/// call and the matching pretty_ua_end_content() -- none of those
/// functions need to change at all. A no-op (content just draws normally,
/// untagged) if the stack is empty or the document was never started --
/// always safe to call defensively around a block even before that
/// block's enclosing function has a pretty_ua_push() of its own yet.
void pretty_ua_begin_content (HPDF_Page page);
void pretty_ua_end_content (HPDF_Page page);

/// Same shape as begin_content()/end_content(), but marks the wrapped
/// block as a PDF/UA Artifact (decorative, outside the structure tree)
/// instead of real content -- for page borders, the migrate logo, and
/// other marks that convey no information of their own.
void pretty_ua_begin_artifact (HPDF_Page page);
void pretty_ua_end_artifact (HPDF_Page page);

/// Called from pdf_new_page() itself (the one place in pretty.c that ever
/// creates a page, whether starting a deliberate new section or handling
/// a plain page-overflow continuation from pdf_advance()/
/// pdf_page_advance_or_not()), right after the new page exists but
/// *before* pdf_new_page()'s own border/header/footer artifact span. A
/// no-op (as it is for the overwhelming majority of pdf_new_page() calls)
/// when nothing is open -- the ordinary case of starting a fresh,
/// independent page.
///
/// If a marked-content span (from begin_content()) or an Artifact span
/// (from begin_artifact()) happens to be open -- a large section's body
/// left one open across what turned out to be a page overflow, rather
/// than closing and reopening it around every individual drawing call --
/// this ends it on `old_page` and remembers what to reopen. The
/// structure element itself is popped, not just suspended:
/// HPDF_UA_BeginMarkedContent() hard-enforces "one struct element's
/// content lives on only one page" and rejects a second page for the
/// same element (an earlier version of this pair tried to reopen the
/// very same element, which the library correctly refused -- and that
/// refusal path triggered pretty.c's own HPDF error handler, which
/// unconditionally longjmp()s out, corrupting this module's state badly
/// enough to crash; caught by this project's own verification, not
/// shipped). Must be paired with a matching pretty_ua_resume_after_new_page()
/// call once the new page's own decoration is done -- see that function.
void pretty_ua_suspend_for_new_page (HPDF_Page old_page);

/// Matching second half of pretty_ua_suspend_for_new_page(): call once
/// `new_page`'s own border/header/footer artifact span (if any) has
/// already been fully opened *and closed* -- reopening any content
/// (as a brand-new sibling element of the same role the suspended one
/// had -- e.g. two Sects, one per page, rather than one spanning both; a
/// real, deliberate trade-off, not a silent bug, and a common real-world
/// PDF/UA pattern besides) or Artifact span that pretty_ua_suspend_for_new_page()
/// had to close, now on `new_page`. A no-op if nothing was suspended.
void pretty_ua_resume_after_new_page (HPDF_Page new_page);

/// Temporarily closes whatever marked-content span is currently open (if
/// any), so a nested structure element -- e.g. a real Table pushed for a
/// matrix embedded partway through a larger Sect's own still-open
/// content span, such as pdf_print_options()'s connection/distance-table
/// calls -- can open its own span without violating "no nesting of open
/// spans" (see pretty_ua_begin_content()). Despite the name (kept for
/// the common case, and to distinguish this pair's call sites from
/// pretty_ua_suspend_for_new_page()'s), this does *not* assume the page
/// stays the same: an earlier version of this pair did, and crashed the
/// first time the nested element's own drawing happened to overflow onto
/// a new page in between -- HPDF_UA_BeginMarkedContent() hard-enforces
/// "one struct element's content lives on only one page", exactly the
/// same rule pretty_ua_suspend_for_new_page() exists to respect, and
/// this pair is now built on the very same mechanism (see pretty_ua.c).
/// Must be paired with a matching pretty_ua_resume_content_same_page()
/// call once the nested element is fully popped back off the stack.
void pretty_ua_suspend_content_same_page (HPDF_Page page);

/// Matching second half of pretty_ua_suspend_content_same_page(). `page`
/// is simply whatever the caller's own `page` variable holds right now
/// -- it does not need to be the same page passed to the suspend call
/// (see that function's comment for why). Reopens whatever was
/// suspended as a brand-new sibling element of the same role -- e.g. two
/// Sects around one Table, rather than one Sect spanning both sides of
/// it -- a real, deliberate trade-off, not a silent bug, and a common
/// real-world PDF/UA pattern besides. A no-op if nothing was suspended.
void pretty_ua_resume_content_same_page (HPDF_Page page);

#endif /*MIGRATE_PRETTY_UA_H*/
