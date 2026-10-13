/* Synced from libharu-pdfua commit 6d0b8eda6670 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
/*
 * hpdf_ua_metadata.c -- Milestone 4: XMP metadata stream (/Metadata).
 *
 * A minimal, purpose-built XMP writer, deliberately NOT a reuse of
 * libharu's own PDF/A HPDF_PDFA_AddXmpMetadata() -- see hpdf_ua.h's
 * comment on HPDF_UA_AddMetadata() for why (it unconditionally
 * (re)creates /MarkInfo and /StructTreeRoot, which would collide with
 * this project's already-tagged, non-empty tree). The XMP packet
 * wrapper/RDF boilerplate below mirrors hpdf_pdfa.c's own HEADER/FOOTER
 * strings (same proven format), narrowed to just dc:title and a
 * PDF/UA-1 identification block (pdfuaid:part) -- the two pieces this
 * project actually needs.
 */
#include <string.h>
#include "hpdf_ua_private.h"
#include "hpdf_utils.h"

#define XMP_HEADER \
    "<?xpacket begin='' id='W5M0MpCehiHzreSzNTczkc9d'?>" \
    "<x:xmpmeta xmlns:x='adobe:ns:meta/'>" \
    "<rdf:RDF xmlns:rdf='http://www.w3.org/1999/02/22-rdf-syntax-ns#'>" \
    "<rdf:Description rdf:about=''" \
    " xmlns:dc='http://purl.org/dc/elements/1.1/'" \
    " xmlns:pdfuaid='http://www.aiim.org/pdfua/ns/id/'>"
#define XMP_TITLE_STARTTAG "<dc:title><rdf:Alt><rdf:li xml:lang=\"x-default\">"
#define XMP_TITLE_ENDTAG   "</rdf:li></rdf:Alt></dc:title>"
#define XMP_PDFUAID_PART1  "<pdfuaid:part>1</pdfuaid:part>"
#define XMP_FOOTER \
    "</rdf:Description></rdf:RDF></x:xmpmeta><?xpacket end='w'?>"

/* Writes `text`, XML-escaped, to `stream`. Returns HPDF_OK, or the first
 * write failure encountered. Truncates silently past a generous fixed
 * buffer -- acceptable for a document title. */
static HPDF_STATUS
write_xml_escaped (HPDF_Stream stream, const char *text)
{
    char buf[2048];
    size_t out = 0;
    const char *p;

    for (p = text; *p && out + 6 < sizeof (buf); p++) {
        switch (*p) {
            case '&':  memcpy (buf + out, "&amp;", 5); out += 5; break;
            case '<':  memcpy (buf + out, "&lt;", 4);  out += 4; break;
            case '>':  memcpy (buf + out, "&gt;", 4);  out += 4; break;
            default:   buf[out++] = *p; break;
        }
    }
    buf[out] = '\0';

    return HPDF_Stream_WriteStr (stream, buf);
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_AddMetadata (HPDF_Doc pdf)
{
    HPDF_Dict xmp;
    const char *title;
    HPDF_STATUS ret = HPDF_OK;

    if (!HPDF_HasDoc (pdf))
        return HPDF_INVALID_DOCUMENT;

    if (HPDF_Dict_GetItem (pdf->catalog, "Metadata", HPDF_OCLASS_DICT))
        return HPDF_OK; /* idempotent */

    xmp = HPDF_DictStream_New (pdf->mmgr, pdf->xref); /* self-registers in the xref */
    if (!xmp)
        return HPDF_CheckError (&pdf->error);

    ret += HPDF_Dict_AddName (xmp, "Type", "Metadata");
    ret += HPDF_Dict_AddName (xmp, "Subtype", "XML");

    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_HEADER);

    title = (const char *) HPDF_GetInfoAttr (pdf, HPDF_INFO_TITLE);
    if (title) {
        ret += HPDF_Stream_WriteStr (xmp->stream, XMP_TITLE_STARTTAG);
        ret += write_xml_escaped (xmp->stream, title);
        ret += HPDF_Stream_WriteStr (xmp->stream, XMP_TITLE_ENDTAG);
    }

    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_PDFUAID_PART1);
    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_FOOTER);

    if (ret != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    if (HPDF_Dict_Add (pdf->catalog, "Metadata", xmp) != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return HPDF_OK;
}
