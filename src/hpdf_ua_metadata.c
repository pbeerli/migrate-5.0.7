/* Synced from libharu-pdfua commit 73af7b6bb214 via tools/sync_libharu_pdfua.sh -- do not hand-edit; edit upstream (https://github.com/pbeerli/libharu-pdfua) and re-run that script. */
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
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "hpdf_ua_private.h"
#include "hpdf_utils.h"
#include "hpdf_encrypt.h" /* HPDF_MD5* -- for the trailer /ID, see
                           * HPDF_UA_AddMetadataWithPDFA()'s own comment */

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

/* Same as XMP_HEADER, plus xmlns:pdfaid -- used only by
 * HPDF_UA_AddMetadataWithPDFA() below. */
#define XMP_HEADER_PDFA \
    "<?xpacket begin='' id='W5M0MpCehiHzreSzNTczkc9d'?>" \
    "<x:xmpmeta xmlns:x='adobe:ns:meta/'>" \
    "<rdf:RDF xmlns:rdf='http://www.w3.org/1999/02/22-rdf-syntax-ns#'>" \
    "<rdf:Description rdf:about=''" \
    " xmlns:dc='http://purl.org/dc/elements/1.1/'" \
    " xmlns:pdfuaid='http://www.aiim.org/pdfua/ns/id/'" \
    " xmlns:pdfaid='http://www.aiim.org/pdfa/ns/id/'>"

/* A real PDF/A Extension Schema declaration (ISO 19005-3 Annex E)
 * describing the `pdfuaid` namespace -- required because `pdfuaid` is
 * not one of PDF/A's own predefined schemas: a strict PDF/A validator
 * otherwise rejects any XMP property in a namespace it doesn't already
 * know about (confirmed directly with veraPDF's own PDF/A-3B flavour,
 * ISO 19005-3:2012 6.6.2.3.1, "properties... shall use either the
 * predefined schemas... or a properly declared extension schema" --
 * before this block was added, HPDF_UA_AddMetadataWithPDFA()'s own
 * pdfuaid:part element failed exactly this check). A second
 * rdf:Description sibling, its own rdf:about='' and namespace
 * declarations, closed and reopened around the main one -- the same
 * multi-Description-per-rdf:RDF shape libharu's own
 * HPDF_PDFA_AddXmpExtension() mechanism uses for its Factur-X example
 * (see vendor/libharu's own demo/pdf_a_conformance.c, upstream, not
 * vendored into this project -- see NOTICE.md). */
#define XMP_PDFUAID_EXTENSION_SCHEMA \
    "<rdf:Description rdf:about=''" \
    " xmlns:pdfaExtension='http://www.aiim.org/pdfa/ns/extension/'" \
    " xmlns:pdfaSchema='http://www.aiim.org/pdfa/ns/schema#'" \
    " xmlns:pdfaProperty='http://www.aiim.org/pdfa/ns/property#'>" \
    "<pdfaExtension:schemas><rdf:Bag><rdf:li rdf:parseType='Resource'>" \
    "<pdfaSchema:schema>PDF/UA identification schema</pdfaSchema:schema>" \
    "<pdfaSchema:namespaceURI>http://www.aiim.org/pdfua/ns/id/</pdfaSchema:namespaceURI>" \
    "<pdfaSchema:prefix>pdfuaid</pdfaSchema:prefix>" \
    "<pdfaSchema:property><rdf:Seq><rdf:li rdf:parseType='Resource'>" \
    "<pdfaProperty:name>part</pdfaProperty:name>" \
    "<pdfaProperty:valueType>Integer</pdfaProperty:valueType>" \
    "<pdfaProperty:category>internal</pdfaProperty:category>" \
    "<pdfaProperty:description>PDF/UA version identifier</pdfaProperty:description>" \
    "</rdf:li></rdf:Seq></pdfaSchema:property>" \
    "</rdf:li></rdf:Bag></pdfaExtension:schemas></rdf:Description>"

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

/* Part digit + conformance letter for each HPDF_PDFA_1A..HPDF_PDFA_3U
 * value, matching hpdf_pdfa.c's own PDFAID_PDFA1A..PDFAID_PDFA3U string
 * constants exactly (same part/conformance pairs, just written as one
 * <pdfaid:part>/<pdfaid:conformance> pair per element instead of a
 * single-tag rdf:Description, so it can sit inside this function's own
 * combined Description alongside pdfuaid). HPDF_PDFA_4/4E/4F use a
 * different, conformance-letter-less shape (see hpdf_pdfa.c's own
 * PDFAID_PDFA4* strings) and are deliberately not supported here yet --
 * no demo in this project needs them, and returning HPDF_OK with no
 * pdfaid element(letting AddMetadataWithPDFA silently downgrade to a
 * plain PDF/UA-1-only declaration for those values) would be a real,
 * silent correctness gap. */
static HPDF_STATUS
write_pdfaid (HPDF_Stream stream, HPDF_PDFAType pdfa_type)
{
    static const struct { HPDF_PDFAType type; char part; char conformance; }
    TABLE[] = {
        { HPDF_PDFA_1A, '1', 'A' }, { HPDF_PDFA_1B, '1', 'B' },
        { HPDF_PDFA_2A, '2', 'A' }, { HPDF_PDFA_2B, '2', 'B' },
        { HPDF_PDFA_2U, '2', 'U' },
        { HPDF_PDFA_3A, '3', 'A' }, { HPDF_PDFA_3B, '3', 'B' },
        { HPDF_PDFA_3U, '3', 'U' },
    };
    size_t i;

    for (i = 0; i < sizeof (TABLE) / sizeof (TABLE[0]); i++) {
        if (TABLE[i].type == pdfa_type) {
            char buf[96];

            (void) snprintf (buf, sizeof (buf),
                    "<pdfaid:part>%c</pdfaid:part>"
                    "<pdfaid:conformance>%c</pdfaid:conformance>",
                    TABLE[i].part, TABLE[i].conformance);
            return HPDF_Stream_WriteStr (stream, buf);
        }
    }

    return HPDF_INVALID_PARAMETER;
}

/* PDF/A (ISO 19005-3:2012 6.1.3) requires the file trailer's /ID entry;
 * PDF/UA-1 does not. libharu has its own equivalent (hpdf_pdfa.c's
 * HPDF_PDFA_GenerateID()), but it is not part of the public API (not
 * declared in hpdf.h) -- this is a real, independent implementation of
 * the same real requirement, not a call to that internal function, and
 * not a modification of vendor/libharu (see NOTICE.md). Does nothing if
 * a trailer /ID already exists (idempotent, matching every other
 * function in this file). */
static HPDF_STATUS
ensure_trailer_id (HPDF_Doc pdf)
{
    HPDF_Array id;
    HPDF_MD5_CTX md5_ctx;
    HPDF_BYTE key[HPDF_MD5_KEY_LEN];
    time_t now;
    char *timestr;

    if (HPDF_Dict_GetItem (pdf->trailer, "ID", HPDF_OCLASS_ARRAY))
        return HPDF_OK;

    id = HPDF_Array_New (pdf->mmgr);
    if (!id)
        return HPDF_CheckError (&pdf->error);
    if (HPDF_Dict_Add (pdf->trailer, "ID", id) != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    now = time (NULL);
    timestr = ctime (&now); /* "Www Mmm dd hh:mm:ss yyyy\n", fixed width */

    HPDF_MD5Init (&md5_ctx);
    HPDF_MD5Update (&md5_ctx, (const HPDF_BYTE *) "libharu-pdfua",
            (HPDF_UINT32) strlen ("libharu-pdfua"));
    HPDF_MD5Update (&md5_ctx, (const HPDF_BYTE *) timestr,
            (HPDF_UINT32) strlen (timestr));
    HPDF_MD5Final (key, &md5_ctx);

    if (HPDF_Array_Add (id, HPDF_Binary_New (pdf->mmgr, key, HPDF_MD5_KEY_LEN)) != HPDF_OK)
        return HPDF_CheckError (&pdf->error);
    if (HPDF_Array_Add (id, HPDF_Binary_New (pdf->mmgr, key, HPDF_MD5_KEY_LEN)) != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return HPDF_OK;
}

HPDF_EXPORT(HPDF_STATUS)
HPDF_UA_AddMetadataWithPDFA (HPDF_Doc pdf, HPDF_PDFAType pdfa_type)
{
    HPDF_Dict xmp;
    const char *title;
    HPDF_STATUS ret = HPDF_OK;

    if (!HPDF_HasDoc (pdf))
        return HPDF_INVALID_DOCUMENT;

    if (HPDF_Dict_GetItem (pdf->catalog, "Metadata", HPDF_OCLASS_DICT))
        return HPDF_OK; /* idempotent */

    xmp = HPDF_DictStream_New (pdf->mmgr, pdf->xref);
    if (!xmp)
        return HPDF_CheckError (&pdf->error);

    ret += HPDF_Dict_AddName (xmp, "Type", "Metadata");
    ret += HPDF_Dict_AddName (xmp, "Subtype", "XML");

    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_HEADER_PDFA);

    title = (const char *) HPDF_GetInfoAttr (pdf, HPDF_INFO_TITLE);
    if (title) {
        ret += HPDF_Stream_WriteStr (xmp->stream, XMP_TITLE_STARTTAG);
        ret += write_xml_escaped (xmp->stream, title);
        ret += HPDF_Stream_WriteStr (xmp->stream, XMP_TITLE_ENDTAG);
    }

    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_PDFUAID_PART1);
    ret += write_pdfaid (xmp->stream, pdfa_type);
    /* Close this Description, add the pdfuaid extension-schema
     * declaration as a sibling Description (see that constant's own
     * comment for why), then close rdf:RDF/x:xmpmeta -- same job
     * XMP_FOOTER does for HPDF_UA_AddMetadata(), with the extension
     * block inserted in the middle. */
    ret += HPDF_Stream_WriteStr (xmp->stream, "</rdf:Description>");
    ret += HPDF_Stream_WriteStr (xmp->stream, XMP_PDFUAID_EXTENSION_SCHEMA);
    ret += HPDF_Stream_WriteStr (xmp->stream,
            "</rdf:RDF></x:xmpmeta><?xpacket end='w'?>");

    if (ret != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    if (HPDF_Dict_Add (pdf->catalog, "Metadata", xmp) != HPDF_OK)
        return HPDF_CheckError (&pdf->error);

    return ensure_trailer_id (pdf);
}
