/* Synced from libharu-pdfua commit b952f3f3a657 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_structure.c -- Milestone 1: structure tree / marked content.
 *
 * Design notes (see hpdf_ua_private.h for the struct layouts these
 * functions build and mutate):
 *
 * - Every structure element's tree linkage (append into its parent's, or
 *   struct_tree_root's, /K array) happens at HPDF_UA_BeginStructureElement()
 *   time, not at End time -- there is deliberately no separate "commit"
 *   step. This means a document remains structurally well-formed even if
 *   an application forgets to call HPDF_UA_EndStructureElement() for some
 *   element (it would just never get freed until HPDF_UA_FreeContext(),
 *   not silently missing from the tree).
 * - A structure element's /K entries are plain MCID integers (not
 *   /MCR marked-content-reference dictionaries), which is valid PDF/UA-1
 *   syntax exactly when every MCID under one element lives on the SAME
 *   page as that element's own /Pg -- which is what
 *   HPDF_UA_BeginMarkedContent() enforces (an element used across more
 *   than one page is rejected, a documented Milestone 1 limitation, not
 *   silently miscompiled).
 * - No /RoleMap entry is written anywhere: every HPDF_UA_StructType maps
 *   directly onto a PDF/UA-1 standard structure type name (see
 *   hpdf_ua_role_name() below), so no custom-to-standard mapping is ever
 *   needed at this milestone -- except HPDF_UA_ROLE_ANNOT ("Annot"), a
 *   real, deliberate exception; see that enum value's own comment in
 *   hpdf_ua.h for why.
 * - Every failure path below routes through HPDF_CheckError()/
 *   HPDF_RaiseError() before returning, even ones this module itself
 *   detects (bad parameters, mismatched Begin/End calls) -- so an
 *   application's registered HPDF_Error_Handler always sees a real
 *   error_no rather than a silent NULL/error status with no diagnostic
 *   trail. Skipping this on early internal returns was a real bug caught
 *   while first testing this file: a NULL came back with nothing printed
 *   anywhere, and the actual cause (see HPDF_UA_NewContext()'s own
 *   version of this issue, if any) took a debugging pass to find instead
 *   of being obvious from the first run's own output.
 */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "hpdf_ua_private.h"
#include "hpdf_utils.h"
#include "hpdf_pages.h"
#include "hpdf_doc.h"

const char *
hpdf_ua_role_name (HPDF_UA_StructType role)
{
    static const char * const NAMES[HPDF_UA_ROLE_EOF] = {
        "Document", "Sect", "Div", "P",
        "H1", "H2", "H3", "H4", "H5", "H6",
        "L", "LI", "Lbl", "LBody",
        "Table", "THead", "TBody", "TFoot", "TR", "TH", "TD",
        "Figure", "Formula", "Caption",
        "Link",
        "Annot",
        "Artifact"
    };

    if (role < 0 || role >= HPDF_UA_ROLE_EOF)
        return NULL;

    return NAMES[role];
}

HPDF_EXPORT(HPDF_UA_StructElem)
HPDF_UA_BeginStructureElement (HPDF_UA_Context ctx, HPDF_UA_StructElem parent,
                                HPDF_UA_StructType role)
{
    HPDF_UA_StructElem elem;
    HPDF_Dict dict;
    HPDF_Array kids;
    const char *role_name;
    void *parent_ref;

    if (!ctx)
        return NULL;

    role_name = hpdf_ua_role_name (role);
    if (!role_name) {
        HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_PARAMETER, 0);
        return NULL;
    }

    dict = HPDF_Dict_New (ctx->pdf->mmgr);
    if (!dict) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }

    if (HPDF_Xref_Add (ctx->pdf->xref, dict) != HPDF_OK) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }

    if (HPDF_Dict_AddName (dict, "Type", "StructElem") != HPDF_OK) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }
    if (HPDF_Dict_AddName (dict, "S", role_name) != HPDF_OK) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }

    /* /P: parent reference. Both a parent element's dict and
     * struct_tree_root are already xref-registered indirect objects
     * (struct_tree_root by HPDF_UA_EnableTagging()), so this is always a
     * real indirect reference, never an accidental inline copy, no
     * matter how many siblings also point /P at the same parent. */
    parent_ref = parent ? (void *) parent->dict : (void *) ctx->struct_tree_root;
    if (HPDF_Dict_Add (dict, "P", parent_ref) != HPDF_OK) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }

    kids = HPDF_Array_New (ctx->pdf->mmgr);
    if (!kids) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }
    if (HPDF_Dict_Add (dict, "K", kids) != HPDF_OK) {
        HPDF_CheckError (&ctx->pdf->error);
        return NULL;
    }

    elem = (HPDF_UA_StructElem) malloc (sizeof (struct _HPDF_UA_StructElem_Rec));
    if (!elem)
        return NULL;

    elem->ctx = ctx;
    elem->dict = dict;
    elem->kids = kids;
    elem->page = NULL;
    elem->next = ctx->elems;
    ctx->elems = elem;

    /* Link into the tree now (see file header comment). */
    if (parent) {
        if (HPDF_Array_Add (parent->kids, dict) != HPDF_OK) {
            HPDF_CheckError (&ctx->pdf->error);
            return NULL;
        }
    } else {
        if (HPDF_Array_Add (ctx->struct_tree_kids, dict) != HPDF_OK) {
            HPDF_CheckError (&ctx->pdf->error);
            return NULL;
        }
    }

    return elem;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_EndStructureElement (HPDF_UA_Context ctx, HPDF_UA_StructElem elem)
{
    /* Deliberately does not free `elem` -- see hpdf_ua_private.h's
     * comment on HPDF_UA_Context_Rec.elems for why: a caller may still
     * legitimately call HPDF_UA_SetAlternateText()/
     * HPDF_UA_SetTableHeaderScope() on it afterward. All elements are
     * freed together by HPDF_UA_FreeContext(). This function exists for
     * call-site symmetry with HPDF_UA_BeginStructureElement() and to
     * validate the handle. */
    if (!ctx || !elem)
        return HPDF_INVALID_PARAMETER;
    if (elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;

    return HPDF_OK;
}

/* Shared body for HPDF_UA_BeginMarkedContent()/_WithActualText() below --
 * identical except for one extra `/ActualText <...>` entry written into
 * the BDC operand dictionary when `actual_text` is non-NULL. Kept as one
 * function (not two independent copies) so the MCID/ParentTree
 * bookkeeping above the BDC write itself can never drift between the
 * two public entry points. */
static HPDF_STATUS
BeginMarkedContentImpl (HPDF_UA_Context ctx, HPDF_Page page,
                         HPDF_UA_StructElem elem, const char *actual_text)
{
    HPDF_UA_PageEntry page_entry;
    HPDF_BOOL first_use_of_elem;
    HPDF_INT32 mcid;
    HPDF_PageAttr attr;

    if (!ctx || !page || !elem || elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;

    if (!HPDF_Page_Validate (page))
        return HPDF_INVALID_PAGE;

    if (elem->page && elem->page != page) {
        /* Documented Milestone 1 limitation: one struct element's content
         * may live on only one page (see this function's header comment
         * in hpdf_ua.h and the /K-as-plain-MCID design note above). */
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_OPERATION, 0);
    }

    page_entry = hpdf_ua_find_or_create_page_entry (ctx, page);
    if (!page_entry)
        return HPDF_CheckError (&ctx->pdf->error);

    if (page_entry->open_elem) {
        /* No nesting in Milestone 1: at most one BDC..EMC span open per
         * page at a time. */
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_OPERATION, 0);
    }

    first_use_of_elem = (HPDF_BOOL) (elem->page == NULL);
    elem->page = page;

    /* /Pg: which page this element's content is on -- pages are always
     * already xref-registered (HPDF_AddPage() does this), so this is
     * always a real indirect reference. Set once, on first use (the
     * check above already guarantees a second use only ever targets the
     * same page, so there is nothing to update on a later call). */
    if (first_use_of_elem) {
        if (HPDF_Dict_Add (elem->dict, "Pg", page) != HPDF_OK)
            return HPDF_CheckError (&ctx->pdf->error);
    }

    mcid = (HPDF_INT32) HPDF_Array_Items (page_entry->mcid_refs);

    /* mcid_refs[mcid] = this element's dict (an indirect reference,
     * since it's xref-registered) -- this is the actual /ParentTree
     * linkage: given a page's /StructParents key and an MCID, a reader
     * resolves ParentTree/Nums[2*key+1][mcid] to find the owning
     * structure element. */
    if (HPDF_Array_Add (page_entry->mcid_refs, elem->dict) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    /* elem's own /K also records this MCID -- see the /K-as-plain-integer
     * design note in this file's header comment. */
    if (HPDF_Array_AddNumber (elem->kids, mcid) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    page_entry->open_elem = elem;

    attr = (HPDF_PageAttr) page->attr;

    /* /<RoleName> << /MCID <n> >> BDC -- the BDC tag conventionally
     * matches the associated structure element's own /S name (best
     * practice; the actual, spec-required linkage is via MCID + Pg +
     * ParentTree, not the tag name itself). Re-fetch /S rather than
     * threading the role enum through this function's own parameters,
     * so the tag written is always exactly what BeginStructureElement()
     * actually recorded. */
    {
        HPDF_Name s_name;
        const char *tag;

        s_name = (HPDF_Name) HPDF_Dict_GetItem (elem->dict, "S",
                HPDF_OCLASS_NAME);
        tag = s_name ? s_name->value : "Span";

        if (HPDF_Stream_WriteEscapeName (attr->stream, tag) != HPDF_OK)
            return HPDF_CheckError (&ctx->pdf->error);
    }

    if (HPDF_Stream_WriteStr (attr->stream, " <</MCID ") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Stream_WriteInt (attr->stream, mcid) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    if (actual_text) {
        /* /ActualText on the marked-content properties dict itself (PDF
         * 32000-1 14.9.4), not just on the structure element -- some
         * checkers/AT walk the content stream's own MC property lists
         * for a text equivalent rather than (or in addition to) a
         * structure element's /ActualText, which is why
         * HPDF_UA_SetActualText() alone (struct-element-only) turned out
         * not to be enough for one real consumer (see that function's
         * own comment history). HPDF_String_Write() gives us the exact
         * same UTF-16BE-with-BOM (or plain PDFDocEncoding, if no "UTF-8"
         * encoder is registered) encoding HPDF_UA_SetActualText() uses,
         * written directly into the content stream instead of into an
         * indirect object -- so this temporary HPDF_String is never
         * added to any dict/xref and must be freed here once written. */
        HPDF_Encoder encoder = HPDF_Doc_FindEncoder (ctx->pdf, "UTF-8");
        HPDF_String s = HPDF_String_New (ctx->pdf->mmgr, actual_text, encoder);

        if (!s)
            return HPDF_CheckError (&ctx->pdf->error);

        if (HPDF_Stream_WriteStr (attr->stream, " /ActualText ") != HPDF_OK) {
            HPDF_String_Free (s);
            return HPDF_CheckError (&ctx->pdf->error);
        }
        if (HPDF_String_Write (s, attr->stream, NULL) != HPDF_OK) {
            HPDF_String_Free (s);
            return HPDF_CheckError (&ctx->pdf->error);
        }
        HPDF_String_Free (s);
    }

    if (HPDF_Stream_WriteStr (attr->stream, ">> BDC\012") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_BeginMarkedContent (HPDF_UA_Context ctx, HPDF_Page page,
                             HPDF_UA_StructElem elem)
{
    return BeginMarkedContentImpl (ctx, page, elem, NULL);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_BeginMarkedContentWithActualText (HPDF_UA_Context ctx, HPDF_Page page,
                                           HPDF_UA_StructElem elem,
                                           const char *actual_text)
{
    if (!actual_text)
        return HPDF_INVALID_PARAMETER;
    return BeginMarkedContentImpl (ctx, page, elem, actual_text);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_EndMarkedContent (HPDF_UA_Context ctx, HPDF_Page page)
{
    HPDF_UA_PageEntry page_entry;
    HPDF_PageAttr attr;

    if (!ctx || !page)
        return HPDF_INVALID_PARAMETER;

    for (page_entry = ctx->pages; page_entry; page_entry = page_entry->next) {
        if (page_entry->page == page)
            break;
    }

    if (!page_entry || !page_entry->open_elem) {
        /* EMC with no matching BDC -- a real caller bug, reported loudly
         * rather than silently emitting an unbalanced operator. */
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_OPERATION, 0);
    }

    attr = (HPDF_PageAttr) page->attr;
    if (HPDF_Stream_WriteStr (attr->stream, "EMC\012") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    page_entry->open_elem = NULL;

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetAlternateText (HPDF_UA_Context ctx, HPDF_UA_StructElem elem,
                           const char *alt_text)
{
    HPDF_String s;

    if (!ctx || !elem || !alt_text || elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;

    s = HPDF_String_New (ctx->pdf->mmgr, alt_text, NULL);
    if (!s)
        return HPDF_CheckError (&ctx->pdf->error);

    HPDF_Dict_RemoveElement (elem->dict, "Alt");

    if (HPDF_Dict_Add (elem->dict, "Alt", s) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetTableHeaderScope (HPDF_UA_Context ctx, HPDF_UA_StructElem th_elem,
                             HPDF_UA_TableScope scope)
{
    static const char * const SCOPE_NAMES[] = { "Row", "Column", "Both" };
    HPDF_Dict attr_dict;
    HPDF_Array attrs;

    if (!ctx || !th_elem || th_elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;
    if (scope != HPDF_UA_SCOPE_ROW && scope != HPDF_UA_SCOPE_COLUMN &&
            scope != HPDF_UA_SCOPE_BOTH)
        return HPDF_INVALID_PARAMETER;

    /* /Scope lives in the element's attribute dictionary under the
     * standard "Table" attribute owner -- PDF 32000-1 14.8.5.4. */
    attrs = (HPDF_Array) HPDF_Dict_GetItem (th_elem->dict, "A", HPDF_OCLASS_ARRAY);
    if (!attrs) {
        attrs = HPDF_Array_New (ctx->pdf->mmgr);
        if (!attrs)
            return HPDF_CheckError (&ctx->pdf->error);
        if (HPDF_Dict_Add (th_elem->dict, "A", attrs) != HPDF_OK)
            return HPDF_CheckError (&ctx->pdf->error);
    }

    attr_dict = HPDF_Dict_New (ctx->pdf->mmgr);
    if (!attr_dict)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Dict_AddName (attr_dict, "O", "Table") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Dict_AddName (attr_dict, "Scope", SCOPE_NAMES[scope]) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    if (HPDF_Array_Add (attrs, attr_dict) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetTableDataHeaders (HPDF_UA_Context ctx, HPDF_UA_StructElem td_elem,
                             HPDF_UA_StructElem *th_elems, HPDF_UINT count)
{
    /* STUB: needs each referenced /TH element to carry a unique /ID, and
     * this /Headers array to reference those IDs by name (PDF 32000-2
     * 14.8.4.5.2) -- more machinery than /Scope alone, deliberately not
     * built out in this milestone's first pass since Migrate's own
     * tables (the motivating consumer) are all simple/regular tables
     * where /Scope is sufficient. Left for a real irregular-table need. */
    (void) ctx;
    (void) td_elem;
    (void) th_elems;
    (void) count;
    return HPDF_UA_NOT_YET_IMPLEMENTED;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_BeginArtifact (HPDF_UA_Context ctx, HPDF_Page page)
{
    return HPDF_UA_BeginArtifactOfType (ctx, page, HPDF_UA_ARTIFACT_NONE);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_BeginArtifactOfType (HPDF_UA_Context ctx, HPDF_Page page,
                              HPDF_UA_ArtifactType type)
{
    static const char * const OPERATORS[] = {
        "/Artifact BMC\012",
        "/Artifact <</Type /Pagination>> BDC\012",
        "/Artifact <</Type /Layout>> BDC\012",
        "/Artifact <</Type /Page>> BDC\012",
        "/Artifact <</Type /Background>> BDC\012"
    };
    HPDF_UA_PageEntry page_entry;
    HPDF_PageAttr attr;

    if (ctx && (type < HPDF_UA_ARTIFACT_NONE
                || type > HPDF_UA_ARTIFACT_BACKGROUND))
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_PARAMETER, 0);

    if (!ctx || !page)
        return HPDF_INVALID_PARAMETER;

    if (!HPDF_Page_Validate (page))
        return HPDF_INVALID_PAGE;

    /* Reuses the same page-entry bookkeeping as marked content (assigns
     * a /StructParents key and an empty mcid_refs array even for a page
     * that turns out to hold only artifacts) -- harmless, and simpler
     * than a second, artifact-only page-tracking path. */
    page_entry = hpdf_ua_find_or_create_page_entry (ctx, page);
    if (!page_entry)
        return HPDF_CheckError (&ctx->pdf->error);

    if (page_entry->open_elem || page_entry->artifact_open) {
        /* No nesting, and no overlap with an open marked-content span --
         * same rule as HPDF_UA_BeginMarkedContent(). */
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_OPERATION, 0);
    }

    page_entry->artifact_open = HPDF_TRUE;

    attr = (HPDF_PageAttr) page->attr;

    /* /Artifact BMC -- BMC (not BDC) is the correct one-operand marked-
     * content operator for "no properties dict" (PDF 32000-1 14.6.2);
     * BDC always requires a second (dict or name) operand. Writing
     * "/Artifact BDC" with no second operand is malformed PDF -- caught
     * empirically, not by inspection: it silently produced a real
     * veraPDF failure (ISO 14289-1:2014 7.1/3, "content shall be marked
     * as Artifact or tagged as real content") on exactly the wrapped
     * region, because a malformed BDC token is not recognized as valid
     * tagging at all. The typed forms carry a property dict, so they do
     * use BDC. */
    if (HPDF_Stream_WriteStr (attr->stream, OPERATORS[type]) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_EndArtifact (HPDF_UA_Context ctx, HPDF_Page page)
{
    HPDF_UA_PageEntry page_entry;
    HPDF_PageAttr attr;

    if (!ctx || !page)
        return HPDF_INVALID_PARAMETER;

    for (page_entry = ctx->pages; page_entry; page_entry = page_entry->next) {
        if (page_entry->page == page)
            break;
    }

    if (!page_entry || !page_entry->artifact_open) {
        return HPDF_RaiseError (&ctx->pdf->error, HPDF_INVALID_OPERATION, 0);
    }

    attr = (HPDF_PageAttr) page->attr;
    if (HPDF_Stream_WriteStr (attr->stream, "EMC\012") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    page_entry->artifact_open = HPDF_FALSE;

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetActualText (HPDF_UA_Context ctx, HPDF_UA_StructElem elem,
                        const char *actual_text)
{
    HPDF_String s;
    HPDF_Encoder encoder;

    if (!ctx || !elem || !actual_text || elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;

    /* /ActualText must be a PDF text string (PDF 32000-1 14.9.4 + 7.9.2.2):
     * PDFDocEncoding bytes, or UTF-16BE with a U+FEFF byte-order mark for
     * anything PDFDocEncoding can't represent -- which is exactly the
     * motivating case here (Greek/math glyphs drawn through a symbol/CID
     * font with no honest PDFDocEncoding spelling). HPDF_String_New()
     * already does that UTF-16BE-with-BOM encoding for us, the same way
     * this library's own Type0/CID text draws do, given a real encoder;
     * reuse whichever one the caller already registered under "UTF-8" via
     * HPDF_UseUTFEncodings() (mirroring exactly how a caller reaches the
     * Identity-H font it drew the glyph through in the first place).
     * HPDF_Doc_FindEncoder() (unlike HPDF_GetEncoder()) just returns NULL
     * with no error side effect when no such encoder is registered, so a
     * caller passing plain ASCII actual_text with no UTF-8 encoder set up
     * falls back to a plain PDFDocEncoding string -- matching
     * HPDF_UA_SetAlternateText()'s own (encoder-less) behavior. */
    encoder = HPDF_Doc_FindEncoder (ctx->pdf, "UTF-8");

    s = HPDF_String_New (ctx->pdf->mmgr, actual_text, encoder);
    if (!s)
        return HPDF_CheckError (&ctx->pdf->error);

    HPDF_Dict_RemoveElement (elem->dict, "ActualText");

    if (HPDF_Dict_Add (elem->dict, "ActualText", s) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_TagAnnotation (HPDF_UA_Context ctx, HPDF_Page page,
                        HPDF_UA_StructElem elem, HPDF_Annotation annot)
{
    HPDF_Dict objr;
    HPDF_INT32 key;

    if (!ctx || !page || !elem || !annot || elem->ctx != ctx)
        return HPDF_INVALID_PARAMETER;

    if (!HPDF_Page_Validate (page))
        return HPDF_INVALID_PAGE;

    /* /OBJR: {Type: /OBJR, Pg: page, Obj: annot} -- the structure-tree
     * kid type used to "contain" an object that isn't page content (an
     * annotation or an XObject), since it has no content-stream position
     * of its own to wrap in BDC/EMC (PDF 32000-1 14.7.4.3). `annot` is
     * always already xref-registered (every HPDF_Page_Create*Annot()
     * constructor does this itself), so /Obj here is always a real
     * indirect reference. */
    objr = HPDF_Dict_New (ctx->pdf->mmgr);
    if (!objr)
        return HPDF_CheckError (&ctx->pdf->error);

    if (HPDF_Dict_AddName (objr, "Type", "OBJR") != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Dict_Add (objr, "Pg", page) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Dict_Add (objr, "Obj", annot) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    if (HPDF_Array_Add (elem->kids, objr) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    /* /StructParent: like a page's /StructParents key, but for a single
     * object rather than an array indexed by MCID (PDF 32000-1
     * 14.7.4.4) -- the /ParentTree /Nums entry for this key is `elem`'s
     * dict directly, not an array of dicts the way a page's entry is.
     * Drawn from the same monotonic counter
     * hpdf_ua_find_or_create_page_entry() uses for pages, so the two
     * key domains never collide even though they share one /Nums tree. */
    key = (HPDF_INT32) (HPDF_Array_Items (ctx->parent_tree_nums) / 2);

    if (HPDF_Array_AddNumber (ctx->parent_tree_nums, key) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);
    if (HPDF_Array_Add (ctx->parent_tree_nums, elem->dict) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    if (HPDF_Dict_AddNumber (annot, "StructParent", key) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    /* /F: Print (bit 3, value 4) set, NoView (bit 6, value 32) clear --
     * ISO 14289-1:2014 7.18 requires this for every annotation included
     * in the logical structure (i.e. every annotation that isn't itself
     * marked as an Artifact). None of libharu's own Link/URI-link
     * annotation constructors set /F at all, so this is normally a
     * fresh add; HPDF_Dict_RemoveElement() first makes it idempotent
     * (safe to call again, or after application code already set some
     * other /F value) rather than failing a second HPDF_Dict_Add(). */
    HPDF_Dict_RemoveElement (annot, "F");
    if (HPDF_Dict_AddNumber (annot, "F", 4) != HPDF_OK)
        return HPDF_CheckError (&ctx->pdf->error);

    /* /Contents: ISO 14289-1:2014 7.18.5 requires link annotations to
     * carry their own alternate description via /Contents (PDF 32000-1
     * 14.9.3) -- a genuinely separate requirement from the /OBJR
     * structure-tree association above, and NOT satisfied by /Alt on
     * `elem` alone (confirmed by a real veraPDF failure on this exact
     * clause the first time this project's own annotation demo was
     * validated). Reuse elem's /Alt text as /Contents when present,
     * rather than requiring a second, separately-worded description --
     * they describe the same link, so one real piece of text serves
     * both purposes. */
    {
        HPDF_String alt_str = (HPDF_String) HPDF_Dict_GetItem (elem->dict,
                "Alt", HPDF_OCLASS_STRING);

        if (alt_str && alt_str->value) {
            HPDF_String contents = HPDF_String_New (ctx->pdf->mmgr,
                    (const char *) alt_str->value, NULL);
            if (!contents)
                return HPDF_CheckError (&ctx->pdf->error);

            HPDF_Dict_RemoveElement (annot, "Contents");
            if (HPDF_Dict_Add (annot, "Contents", contents) != HPDF_OK)
                return HPDF_CheckError (&ctx->pdf->error);
        }
    }

    return HPDF_OK;
}
