package header

// Reserved keyword constants from FITS Standard v4.0, §4.4 and Appendix C.
//
// These names are provided so callers can compare against a single source of
// truth rather than sprinkling string literals. The list is not exhaustive;
// we add as needed.

// Mandatory structural keywords.
const (
	KeySimple   = "SIMPLE"
	KeyBitpix   = "BITPIX"
	KeyNaxis    = "NAXIS"
	KeyXtension = "XTENSION"
	KeyEnd      = "END"
	KeyExtname  = "EXTNAME"
	KeyExtver   = "EXTVER"
	KeyExtlevel = "EXTLEVEL"
	KeyPcount   = "PCOUNT"
	KeyGcount   = "GCOUNT"
	KeyTfields  = "TFIELDS"
)

// Commentary keywords.
const (
	KeyComment = "COMMENT"
	KeyHistory = "HISTORY"
	KeyBlank   = "" // blank keyword
)

// Scaling keywords.
const (
	KeyBscale = "BSCALE"
	KeyBzero  = "BZERO"
	KeyBlankV = "BLANK" // "BLANK" keyword (distinct from blank keyword constant)
	KeyBunit  = "BUNIT"
)

// XTENSION values.
const (
	XtensionImage    = "IMAGE"
	XtensionBinTable = "BINTABLE"
	XtensionTable    = "TABLE"
)

// Checksum keywords.
const (
	KeyChecksum = "CHECKSUM"
	KeyDatasum  = "DATASUM"
)
