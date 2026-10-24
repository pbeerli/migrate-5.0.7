/* Synced from libharu-pdfua commit 73af7b6bb214 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_context.c -- HPDF_UA_Context lifecycle (Milestone 1).
 *
 * See hpdf_ua.h's comment above the Context section for why this exists
 * (persistent tagging bookkeeping needs an explicit owner; this project
 * does not modify vendored libharu headers to add fields to HPDF_Doc
 * itself, and a module-global registry keyed by HPDF_Doc would be unsafe
 * across multiple simultaneously-open documents and doc lifetimes).
 *
 * This project's own bookkeeping structs (HPDF_UA_Context_Rec,
 * HPDF_UA_PageEntry_Rec) are plain C structs allocated with the standard
 * library, not libharu's HPDF_MMgr -- they are not PDF objects and are
 * never serialized; only the real PDF dict/array objects they point to
 * (created via HPDF_Dict_New()/HPDF_Array_New() against pdf->mmgr, exactly
 * like the rest of this module) are.
 */
#include <stdlib.h>
#include "hpdf_ua_private.h"
#include "hpdf_utils.h"

HPDF_EXPORT(HPDF_UA_Context)
HPDF_UA_NewContext (HPDF_Doc pdf)
{
    HPDF_UA_Context ctx;
    HPDF_Dict struct_tree_root;
    HPDF_Array struct_tree_kids;
    HPDF_Dict parent_tree;
    HPDF_Array parent_tree_nums;

    if (!HPDF_HasDoc (pdf))
        return NULL;

    if (HPDF_UA_EnableTagging (pdf) != HPDF_OK)
        return NULL;

    struct_tree_root = (HPDF_Dict) HPDF_Dict_GetItem (pdf->catalog,
            "StructTreeRoot", HPDF_OCLASS_DICT);
    if (!struct_tree_root)
        return NULL; /* HPDF_UA_EnableTagging() should have created this */

    struct_tree_kids = (HPDF_Array) HPDF_Dict_GetItem (struct_tree_root, "K",
            HPDF_OCLASS_ARRAY);
    if (!struct_tree_kids)
        return NULL;

    /* /ParentTree: created once, reused if a context was already set up
     * for this document (re-fetch rather than assume this is the first
     * HPDF_UA_NewContext() call for this pdf). */
    parent_tree = (HPDF_Dict) HPDF_Dict_GetItem (struct_tree_root,
            "ParentTree", HPDF_OCLASS_DICT);
    if (!parent_tree) {
        parent_tree = HPDF_Dict_New (pdf->mmgr);
        if (!parent_tree)
            return NULL;
        if (HPDF_Xref_Add (pdf->xref, parent_tree) != HPDF_OK)
            return NULL;

        parent_tree_nums = HPDF_Array_New (pdf->mmgr);
        if (!parent_tree_nums)
            return NULL;
        if (HPDF_Dict_Add (parent_tree, "Nums", parent_tree_nums) != HPDF_OK)
            return NULL;

        if (HPDF_Dict_Add (struct_tree_root, "ParentTree", parent_tree)
                != HPDF_OK)
            return NULL;
    } else {
        parent_tree_nums = (HPDF_Array) HPDF_Dict_GetItem (parent_tree,
                "Nums", HPDF_OCLASS_ARRAY);
        if (!parent_tree_nums)
            return NULL;
    }

    ctx = (HPDF_UA_Context) malloc (sizeof (struct _HPDF_UA_Context_Rec));
    if (!ctx)
        return NULL;

    ctx->pdf = pdf;
    ctx->struct_tree_root = struct_tree_root;
    ctx->struct_tree_kids = struct_tree_kids;
    ctx->parent_tree = parent_tree;
    ctx->parent_tree_nums = parent_tree_nums;
    ctx->pages = NULL;
    ctx->elems = NULL;

    return ctx;
}

HPDF_EXPORT(void)
HPDF_UA_FreeContext (HPDF_UA_Context ctx)
{
    HPDF_UA_PageEntry page_entry, next_page_entry;
    HPDF_UA_StructElem elem, next_elem;

    if (!ctx)
        return;

    /* Frees only this project's own bookkeeping list nodes -- the PDF
     * dict/array objects they point to (mcid_refs, struct elem dicts,
     * struct_tree_root, parent_tree) remain owned by pdf->mmgr/pdf->xref
     * exactly as any other libharu object, and are torn down normally by
     * HPDF_Free(pdf), not by this call. */
    page_entry = ctx->pages;
    while (page_entry) {
        next_page_entry = page_entry->next;
        free (page_entry);
        page_entry = next_page_entry;
    }

    elem = ctx->elems;
    while (elem) {
        next_elem = elem->next;
        free (elem);
        elem = next_elem;
    }

    free (ctx);
}

HPDF_UA_PageEntry
hpdf_ua_find_or_create_page_entry (HPDF_UA_Context ctx, HPDF_Page page)
{
    HPDF_UA_PageEntry entry;
    HPDF_Array mcid_refs;
    HPDF_INT32 key;

    for (entry = ctx->pages; entry; entry = entry->next) {
        if (entry->page == page)
            return entry;
    }

    /* First time this page has been tagged: assign it the next
     * /StructParents key (derived from how many [key, arr] pairs already
     * exist in /Nums, so this stays correct even if HPDF_UA_NewContext()
     * re-attached to a /ParentTree a previous context had already
     * started populating). */
    key = (HPDF_INT32) (HPDF_Array_Items (ctx->parent_tree_nums) / 2);

    if (HPDF_Dict_AddNumber (page, "StructParents", key) != HPDF_OK)
        return NULL;

    /* /Tabs /S (Milestone 4): tab/reading order on this page follows
     * structure order -- a real PDF/UA-1 requirement, set automatically
     * here rather than left to each caller to remember for every page it
     * tags. */
    if (HPDF_Dict_AddName (page, "Tabs", "S") != HPDF_OK)
        return NULL;

    mcid_refs = HPDF_Array_New (ctx->pdf->mmgr);
    if (!mcid_refs)
        return NULL;

    /* Eagerly link [key, mcid_refs] into /Nums now, even though
     * mcid_refs is still empty -- later HPDF_Array_Add() calls onto this
     * same array (as more marked-content spans are opened on this page)
     * mutate it in place, and are automatically reflected at
     * HPDF_SaveToFile() time since serialization happens once, at the
     * very end. No separate finalize step is needed. */
    if (HPDF_Array_AddNumber (ctx->parent_tree_nums, key) != HPDF_OK)
        return NULL;
    if (HPDF_Array_Add (ctx->parent_tree_nums, mcid_refs) != HPDF_OK)
        return NULL;

    entry = (HPDF_UA_PageEntry) malloc (sizeof (HPDF_UA_PageEntry_Rec));
    if (!entry)
        return NULL;

    entry->page = page;
    entry->struct_parents_key = key;
    entry->mcid_refs = mcid_refs;
    entry->open_elem = NULL;
    entry->artifact_open = HPDF_FALSE;
    entry->next = ctx->pages;
    ctx->pages = entry;

    return entry;
}
