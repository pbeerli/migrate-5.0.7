/* Synced from libharu-pdfua commit b952f3f3a657 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_private.h -- internal struct definitions shared by this
 * module's own .c files. Not installed/exported; application code should
 * only ever see the opaque handles declared in hpdf_ua.h.
 */
#ifndef _HPDF_UA_PRIVATE_H
#define _HPDF_UA_PRIVATE_H

#include "hpdf_ua.h"

/* One entry per page this context has tagged content on. Owned by, and
 * freed together with, the HPDF_UA_Context that created it -- see
 * hpdf_ua_context.c. Never persists beyond its owning context's lifetime,
 * so there is no cross-document-lifetime stale-pointer risk the way a
 * module-global registry keyed by HPDF_Page would have. */
typedef struct _HPDF_UA_PageEntry_Rec {
    HPDF_Page                       page;
    HPDF_INT32                      struct_parents_key;
    HPDF_Array                      mcid_refs;  /* index = MCID, value = the
                                                  * struct elem dict that owns
                                                  * it (an indirect reference,
                                                  * since the elem dict is
                                                  * xref-registered) */
    HPDF_UA_StructElem               open_elem;  /* non-NULL while a BDC..EMC
                                                  * marked-content span is
                                                  * open on this page; at
                                                  * most one at a time
                                                  * (no nesting) */
    HPDF_BOOL                       artifact_open; /* an /Artifact BDC..EMC
                                                  * span (HPDF_UA_BeginArtifact())
                                                  * is open on this page --
                                                  * mutually exclusive with
                                                  * open_elem, same no-nesting
                                                  * rule */
    struct _HPDF_UA_PageEntry_Rec  *next;
} HPDF_UA_PageEntry_Rec, *HPDF_UA_PageEntry;

struct _HPDF_UA_Context_Rec {
    HPDF_Doc            pdf;
    HPDF_Dict           struct_tree_root;
    HPDF_Array          struct_tree_kids;   /* struct_tree_root's own /K */
    HPDF_Dict           parent_tree;
    HPDF_Array          parent_tree_nums;   /* flat /Nums: [key0 arr0 key1 arr1 ...] --
                                              * a single flat number tree, not a
                                              * balanced /Kids tree; fine for the
                                              * page counts this project's actual
                                              * consumer (Migrate reports) produces,
                                              * documented as a limitation. */
    HPDF_UA_PageEntry   pages;              /* linked list, linear scan --
                                              * page counts are small (see above) */
    HPDF_UA_StructElem  elems;              /* linked list of every struct-elem
                                              * wrapper this context has handed
                                              * out (see struct below); freed in
                                              * bulk by HPDF_UA_FreeContext(),
                                              * deliberately NOT freed by
                                              * HPDF_UA_EndStructureElement() --
                                              * a caller may legitimately still
                                              * want to call
                                              * HPDF_UA_SetAlternateText()/
                                              * HPDF_UA_SetTableHeaderScope() on
                                              * an element after "ending" it,
                                              * and freeing early would make
                                              * that a use-after-free. */
};

/* One node in the logical structure tree. `dict` is xref-registered
 * (indirect) at creation, since every structure element is referenced
 * from at least two places once tagged content exists under it: its
 * parent's /K array, and (once it owns marked content) a page's
 * mcid_refs array -- an object embedded inline could not be shared that
 * way; see docs/pdf_ua_requirements.md and the design notes in
 * hpdf_ua_structure.c. */
struct _HPDF_UA_StructElem_Rec {
    HPDF_UA_Context     ctx;
    HPDF_Dict           dict;
    HPDF_Array          kids;   /* this element's own /K array */
    HPDF_Page           page;   /* the one page this elem's content lives on,
                                  * once first used with
                                  * HPDF_UA_BeginMarkedContent(); NULL until
                                  * then. See that function's docs on the
                                  * single-page-per-element limitation. */
    HPDF_UA_StructElem  next;   /* ctx->elems linked list */
};

/* Internal helper, defined in hpdf_ua_structure.c, used by both files:
 * maps a role enum to its standard PDF structure-type name (also used
 * verbatim as the BDC tag name for marked content under that element,
 * per Tagged-PDF convention). */
const char *
hpdf_ua_role_name (HPDF_UA_StructType role);

/* Internal helper, defined in hpdf_ua_context.c: finds this context's
 * existing page-entry for `page`, or creates one (assigning the next
 * /StructParents key, writing it onto the page dict, creating an empty
 * mcid_refs array, and eagerly linking [key, mcid_refs] into
 * /ParentTree's /Nums so later growth of mcid_refs is automatically
 * reflected at save time with no separate finalize step). Returns NULL
 * only on allocation failure. */
HPDF_UA_PageEntry
hpdf_ua_find_or_create_page_entry (HPDF_UA_Context ctx, HPDF_Page page);

#endif /* _HPDF_UA_PRIVATE_H */
