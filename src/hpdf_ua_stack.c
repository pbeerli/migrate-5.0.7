/* Synced from libharu-pdfua commit 73af7b6bb214 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_stack.c -- the tag stack (0.9.0), see hpdf_ua.h.
 *
 * Generalized from Migrate-n's pretty_ua.c, which wraps its whole tagged
 * PDF report (tables and figures that routinely run over many pages) in
 * exactly this push/pop/suspend/resume discipline. Differences from that
 * original: the state lives in an HPDF_UA_Stack the caller owns instead
 * of module globals, the stack grows without a fixed depth limit, the
 * stack remembers the page its open span is on, and misuse is reported
 * through HPDF_RaiseError() like the rest of this module instead of
 * being ignored.
 *
 * Why suspend pops the element instead of only closing its span:
 * HPDF_UA_BeginMarkedContent() allows one structure element's content on
 * one page only (see hpdf_ua_structure.c's header), and at suspend time
 * there is no way to know whether the gap until resume will cross a page.
 * Always popping, and on resume always pushing a new sibling of the same
 * role, is correct either way, however many pages go by in between.
 *
 * Why the two suspend/resume pairs keep separate state: see the comment
 * above HPDF_UA_Stack_SuspendContent() in hpdf_ua.h. Migrate-n first used
 * one shared slot and lost the tags of every nested table that straddled
 * a page break (found by a veraPDF before/after comparison).
 */
#include <stdlib.h>
#include "hpdf_ua_private.h"
#include "hpdf_utils.h"

typedef struct {
    HPDF_BOOL            content;      /* a content span was suspended */
    HPDF_UA_StructType   content_role; /* ...for an element of this role */
    HPDF_BOOL            artifact;     /* an artifact span was suspended */
    HPDF_UA_ArtifactType artifact_type;
} HPDF_UA_Suspended;

struct _HPDF_UA_Stack_Rec {
    HPDF_UA_Context      ctx;
    HPDF_UA_StructElem  *elems;
    HPDF_UA_StructType  *roles;
    HPDF_UINT            depth;
    HPDF_UINT            capacity;
    HPDF_Page            content_page;  /* non-NULL while a content span
                                          * is open for the top element */
    HPDF_Page            artifact_page; /* non-NULL while an artifact span
                                          * is open */
    HPDF_UA_ArtifactType artifact_type;
    HPDF_UA_Suspended    page_break;    /* SuspendForNewPage/ResumeAfterNewPage */
    HPDF_UA_Suspended    nested;        /* SuspendContent/ResumeContent */
};

static HPDF_STATUS
stack_error (HPDF_UA_Stack stack, HPDF_STATUS code)
{
    return HPDF_RaiseError (&stack->ctx->pdf->error, code, 0);
}

HPDF_EXPORT(HPDF_UA_Stack)
HPDF_UA_NewStack (HPDF_UA_Context ctx)
{
    HPDF_UA_Stack stack;

    if (!ctx)
        return NULL;

    stack = (HPDF_UA_Stack) calloc (1, sizeof (*stack));
    if (!stack) {
        HPDF_RaiseError (&ctx->pdf->error, HPDF_FAILD_TO_ALLOC_MEM, 0);
        return NULL;
    }
    stack->ctx = ctx;
    return stack;
}

HPDF_EXPORT(void)
HPDF_UA_FreeStack (HPDF_UA_Stack stack)
{
    if (!stack)
        return;

    if (stack->artifact_page)
        HPDF_UA_Stack_EndArtifact (stack);
    while (stack->depth > 0)
        HPDF_UA_Stack_Pop (stack);   /* also closes an open content span */

    free (stack->elems);
    free (stack->roles);
    free (stack);
}

HPDF_EXPORT(HPDF_UA_StructElem)
HPDF_UA_Stack_Push (HPDF_UA_Stack stack, HPDF_UA_StructType role)
{
    HPDF_UA_StructElem parent;
    HPDF_UA_StructElem elem;

    if (!stack)
        return NULL;

    if (stack->depth == stack->capacity) {
        HPDF_UINT capacity = stack->capacity ? 2 * stack->capacity : 16;
        HPDF_UA_StructElem *elems = (HPDF_UA_StructElem *)
                realloc (stack->elems, capacity * sizeof (*elems));
        HPDF_UA_StructType *roles;

        if (!elems) {
            stack_error (stack, HPDF_FAILD_TO_ALLOC_MEM);
            return NULL;
        }
        stack->elems = elems;
        roles = (HPDF_UA_StructType *)
                realloc (stack->roles, capacity * sizeof (*roles));
        if (!roles) {
            stack_error (stack, HPDF_FAILD_TO_ALLOC_MEM);
            return NULL;
        }
        stack->roles = roles;
        stack->capacity = capacity;
    }

    parent = stack->depth > 0 ? stack->elems[stack->depth - 1] : NULL;
    elem = HPDF_UA_BeginStructureElement (stack->ctx, parent, role);
    if (!elem)
        return NULL;

    stack->elems[stack->depth] = elem;
    stack->roles[stack->depth] = role;
    stack->depth++;
    return elem;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_Pop (HPDF_UA_Stack stack)
{
    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (stack->depth == 0)
        return stack_error (stack, HPDF_INVALID_OPERATION);

    if (stack->content_page)
        HPDF_UA_Stack_EndContent (stack);

    stack->depth--;
    return HPDF_UA_EndStructureElement (stack->ctx,
            stack->elems[stack->depth]);
}

HPDF_EXPORT(HPDF_UA_StructElem)
HPDF_UA_Stack_Top (HPDF_UA_Stack stack)
{
    if (!stack || stack->depth == 0)
        return NULL;
    return stack->elems[stack->depth - 1];
}

HPDF_EXPORT(HPDF_UINT)
HPDF_UA_Stack_Depth (HPDF_UA_Stack stack)
{
    return stack ? stack->depth : 0;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_SetAlternateText (HPDF_UA_Stack stack, const char *alt_text)
{
    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (stack->depth == 0)
        return stack_error (stack, HPDF_INVALID_OPERATION);
    return HPDF_UA_SetAlternateText (stack->ctx,
            stack->elems[stack->depth - 1], alt_text);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_SetActualText (HPDF_UA_Stack stack, const char *actual_text)
{
    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (stack->depth == 0)
        return stack_error (stack, HPDF_INVALID_OPERATION);
    return HPDF_UA_SetActualText (stack->ctx,
            stack->elems[stack->depth - 1], actual_text);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_SetTableHeaderScope (HPDF_UA_Stack stack,
                                    HPDF_UA_TableScope scope)
{
    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (stack->depth == 0)
        return stack_error (stack, HPDF_INVALID_OPERATION);
    return HPDF_UA_SetTableHeaderScope (stack->ctx,
            stack->elems[stack->depth - 1], scope);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_BeginContentWithActualText (HPDF_UA_Stack stack,
                                           HPDF_Page page,
                                           const char *actual_text)
{
    HPDF_UA_StructElem top;
    HPDF_STATUS status;

    if (!stack || !page)
        return HPDF_INVALID_PARAMETER;
    if (stack->depth == 0 || stack->content_page || stack->artifact_page)
        return stack_error (stack, HPDF_INVALID_OPERATION);

    top = stack->elems[stack->depth - 1];
    if (actual_text)
        status = HPDF_UA_BeginMarkedContentWithActualText (stack->ctx, page,
                top, actual_text);
    else
        status = HPDF_UA_BeginMarkedContent (stack->ctx, page, top);
    if (status == HPDF_OK)
        stack->content_page = page;
    return status;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_BeginContent (HPDF_UA_Stack stack, HPDF_Page page)
{
    return HPDF_UA_Stack_BeginContentWithActualText (stack, page, NULL);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_EndContent (HPDF_UA_Stack stack)
{
    HPDF_Page page;

    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (!stack->content_page)
        return stack_error (stack, HPDF_INVALID_OPERATION);

    page = stack->content_page;
    stack->content_page = NULL;
    return HPDF_UA_EndMarkedContent (stack->ctx, page);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_BeginArtifact (HPDF_UA_Stack stack, HPDF_Page page,
                              HPDF_UA_ArtifactType type)
{
    HPDF_STATUS status;

    if (!stack || !page)
        return HPDF_INVALID_PARAMETER;
    if (stack->artifact_page || stack->content_page)
        return stack_error (stack, HPDF_INVALID_OPERATION);

    status = HPDF_UA_BeginArtifactOfType (stack->ctx, page, type);
    if (status == HPDF_OK) {
        stack->artifact_page = page;
        stack->artifact_type = type;
    }
    return status;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_EndArtifact (HPDF_UA_Stack stack)
{
    HPDF_Page page;

    if (!stack)
        return HPDF_INVALID_PARAMETER;
    if (!stack->artifact_page)
        return stack_error (stack, HPDF_INVALID_OPERATION);

    page = stack->artifact_page;
    stack->artifact_page = NULL;
    return HPDF_UA_EndArtifact (stack->ctx, page);
}

/* Shared engine behind both public suspend/resume pairs, parameterized
 * on which pair's saved state to use. */
static HPDF_STATUS
suspend_and_pop (HPDF_UA_Stack stack, HPDF_UA_Suspended *saved)
{
    HPDF_STATUS status = HPDF_OK;

    if (!stack)
        return HPDF_INVALID_PARAMETER;

    saved->content = (HPDF_BOOL) (stack->content_page != NULL);
    saved->artifact = (HPDF_BOOL) (stack->artifact_page != NULL);

    if (saved->content) {
        saved->content_role = stack->roles[stack->depth - 1];
        status = HPDF_UA_Stack_Pop (stack);   /* ends the span first */
    }
    if (saved->artifact) {
        saved->artifact_type = stack->artifact_type;
        if (status == HPDF_OK)
            status = HPDF_UA_Stack_EndArtifact (stack);
    }
    return status;
}

static HPDF_STATUS
resume_as_new_sibling (HPDF_UA_Stack stack, HPDF_UA_Suspended *saved,
                       HPDF_Page page)
{
    HPDF_STATUS status = HPDF_OK;

    if (!stack || !page)
        return HPDF_INVALID_PARAMETER;

    if (saved->content) {
        saved->content = HPDF_FALSE;
        if (!HPDF_UA_Stack_Push (stack, saved->content_role))
            return HPDF_CheckError (&stack->ctx->pdf->error);
        status = HPDF_UA_Stack_BeginContent (stack, page);
    }
    if (saved->artifact) {
        saved->artifact = HPDF_FALSE;
        if (status == HPDF_OK)
            status = HPDF_UA_Stack_BeginArtifact (stack, page,
                    saved->artifact_type);
    }
    return status;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_SuspendForNewPage (HPDF_UA_Stack stack)
{
    return stack ? suspend_and_pop (stack, &stack->page_break)
                 : HPDF_INVALID_PARAMETER;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_ResumeAfterNewPage (HPDF_UA_Stack stack, HPDF_Page new_page)
{
    return stack ? resume_as_new_sibling (stack, &stack->page_break, new_page)
                 : HPDF_INVALID_PARAMETER;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_SuspendContent (HPDF_UA_Stack stack)
{
    return stack ? suspend_and_pop (stack, &stack->nested)
                 : HPDF_INVALID_PARAMETER;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_Stack_ResumeContent (HPDF_UA_Stack stack, HPDF_Page page)
{
    return stack ? resume_as_new_sibling (stack, &stack->nested, page)
                 : HPDF_INVALID_PARAMETER;
}
