/* Synced from libharu-pdfua commit 6d0b8eda6670 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_docmeta.c -- real, working Milestone-0 additions: document
 * language, DisplayDocTitle, and idempotent MarkInfo/StructTreeRoot setup.
 *
 * This file is this project's own addition (see NOTICE.md), built against
 * libharu's public API plus its internal object headers (hpdf_doc.h,
 * hpdf_objects.h) -- the same headers libharu's own hpdf_pdfa.c uses to do
 * the equivalent kind of catalog-level dictionary manipulation. No vendored
 * libharu source file is modified to make this work.
 */
#include "hpdf_ua.h"
#include "hpdf_utils.h"

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetDocumentLanguage (HPDF_Doc pdf, const char *lang)
{
    HPDF_STATUS ret;
    HPDF_String lang_str;

    if (!HPDF_HasDoc (pdf))
        return HPDF_INVALID_DOCUMENT;

    if (!lang)
        return HPDF_RaiseError (&pdf->error, HPDF_INVALID_PARAMETER, 0);

    lang_str = HPDF_String_New (pdf->mmgr, lang, NULL);
    if (!lang_str)
        return HPDF_CheckError (&pdf->error);

    /* Replace any existing /Lang entry rather than erroring on a second
     * call -- callers may legitimately want to change the language after
     * the fact (e.g. multi-locale report generation). */
    HPDF_Dict_RemoveElement (pdf->catalog, "Lang");

    ret = HPDF_Dict_Add (pdf->catalog, "Lang", lang_str);
    if (ret != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_SetDisplayDocTitle (HPDF_Doc pdf, HPDF_BOOL value)
{
    HPDF_STATUS ret;
    HPDF_Dict preferences;

    if (!HPDF_HasDoc (pdf))
        return HPDF_INVALID_DOCUMENT;

    preferences = (HPDF_Dict) HPDF_Dict_GetItem (pdf->catalog,
            "ViewerPreferences", HPDF_OCLASS_DICT);

    if (!preferences) {
        preferences = HPDF_Dict_New (pdf->mmgr);
        if (!preferences)
            return HPDF_CheckError (&pdf->error);

        ret = HPDF_Dict_Add (pdf->catalog, "ViewerPreferences", preferences);
        if (ret != HPDF_OK)
            return HPDF_CheckError (&pdf->error);
    }

    ret = HPDF_Dict_AddBoolean (preferences, "DisplayDocTitle", value);
    if (ret != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_EnableTagging (HPDF_Doc pdf)
{
    HPDF_STATUS ret = HPDF_OK;
    HPDF_Dict markinfo;
    HPDF_Dict struct_tree_root;
    HPDF_Array k;

    if (!HPDF_HasDoc (pdf))
        return HPDF_INVALID_DOCUMENT;

    /* Idempotent: if a previous call (or the legacy PDF/A path in
     * hpdf_pdfa.c, if the application also uses that) already created
     * these, leave them alone rather than clobbering real content a later
     * milestone may have already attached. */
    markinfo = (HPDF_Dict) HPDF_Dict_GetItem (pdf->catalog, "MarkInfo",
            HPDF_OCLASS_DICT);
    if (!markinfo) {
        markinfo = HPDF_Dict_New (pdf->mmgr);
        if (!markinfo)
            return HPDF_CheckError (&pdf->error);

        ret += HPDF_Dict_Add (pdf->catalog, "MarkInfo", markinfo);
        ret += HPDF_Dict_AddBoolean (markinfo, "Marked", HPDF_TRUE);
    }

    struct_tree_root = (HPDF_Dict) HPDF_Dict_GetItem (pdf->catalog,
            "StructTreeRoot", HPDF_OCLASS_DICT);
    if (!struct_tree_root) {
        struct_tree_root = HPDF_Dict_New (pdf->mmgr);
        if (!struct_tree_root)
            return HPDF_CheckError (&pdf->error);

        /* Must be xref-registered (indirect) BEFORE being added anywhere,
         * unlike the original PDF/A path this was adapted from (which
         * never needed to reference struct_tree_root a second time).
         * Milestone 1's structure elements all point their own /P back at
         * struct_tree_root for top-level elements -- and libharu's object
         * model marks any object HPDF_OTYPE_DIRECT ("owned by exactly one
         * container") the moment it is added anywhere via HPDF_Dict_Add()/
         * HPDF_Array_Add() unless it is already xref-registered first,
         * silently refusing any second HPDF_Dict_Add()/HPDF_Array_Add()
         * call with HPDF_INVALID_OBJECT. Found the hard way, empirically,
         * while first testing Milestone 1: every single
         * HPDF_UA_BeginStructureElement() call for a top-level element
         * failed until this fix, since each one tries to add
         * struct_tree_root as its own /P value. */
        if (HPDF_Xref_Add (pdf->xref, struct_tree_root) != HPDF_OK)
            return HPDF_CheckError (&pdf->error);

        ret += HPDF_Dict_Add (pdf->catalog, "StructTreeRoot", struct_tree_root);
        ret += HPDF_Dict_AddName (struct_tree_root, "Type", "StructTreeRoot");

        /* Milestone 0: still an empty /K array, exactly like the existing
         * PDF/A path -- Milestone 1 is what actually populates this with
         * real structure elements via HPDF_UA_BeginStructureElement(). */
        k = HPDF_Array_New (pdf->mmgr);
        if (!k)
            return HPDF_CheckError (&pdf->error);
        ret += HPDF_Dict_Add (struct_tree_root, "K", k);
    }

    if (ret != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return HPDF_OK;
}
