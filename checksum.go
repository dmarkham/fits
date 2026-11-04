package fits

import (
	"errors"
	"fmt"
	"math"
	"strconv"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/block"
	"github.com/dmarkham/fits/internal/checksum"
)

// The FITS checksum convention (FITS Standard v4.0 Appendix J): DATASUM is
// the 32-bit ones'-complement sum of the data unit (including its block
// padding) written as a decimal string; CHECKSUM is a 16-character ASCII
// encoding chosen so that the ones'-complement sum of the whole HDU
// (header blocks, with the CHECKSUM card as stored, plus the data unit) is
// all ones (negative zero). The data sum and header sum are combined with
// the same end-around-carry addition, as cfitsio's ffpcks/ffvcks do.

// ChecksumStatus verifies an HDU's DATASUM and CHECKSUM keywords against
// the bytes on disk and reports each the way cfitsio's fits_verify_chksum
// does: 1 when the keyword is present and correct, 0 when it is absent (or
// blank), -1 when it is present and wrong. It works for every HDU kind,
// including tile-compressed images (whose checksum covers the binary table
// as stored). The returned error is for I/O problems, and
// ErrUndefinedChecksum for a DATASUM or CHECKSUM card with no value (where
// fits_verify_chksum fails too).
func ChecksumStatus(h HDU) (datasum, hdusum int, err error) {
	rec := hduRecordOf(h)
	if rec == nil || rec.file == nil || rec.file.br == nil {
		return 0, 0, errors.New("fits: ChecksumStatus: HDU is not backed by a readable file")
	}
	hdr, dataSum, err := rec.onDiskSums()
	if err != nil {
		return 0, 0, err
	}
	// The keywords are found and parsed as cfitsio's ffvcks does. It reads
	// CHECKSUM first and then DATASUM, each with ffgcrd, which searches
	// from the card after the last one it read to END and then wraps to
	// the top: CHECKSUM is the first card of that name, and DATASUM the
	// first one after that CHECKSUM card (or, if there is none after it
	// or no CHECKSUM, the first from the top). Names are compared in
	// upper case; the value is unquoted (ffgkys) and empty means absent.
	// DATASUM is read with atof (the longest numeric prefix, so "123abc"
	// is 123 and hex floats are accepted) and cast to unsigned long, which
	// truncates toward zero: it matches when the truncated value equals
	// the data sum (so "123.9" matches 123 and "-0.5", which truncates to
	// 0, matches a zero sum). A value at or below -1, or of 2^32 or more,
	// never matches: cfitsio's cast of a negative value is undefined in C
	// (on x86-64 it yields 2^64 minus the magnitude, which cannot equal a
	// 32-bit sum) and a larger value cannot equal it either (unsigned long
	// is 64 bits). A card with no value at all (not even '') makes
	// cfitsio's ffgkys fail with VALUE_UNDEFINED, and fits_verify_chksum
	// with it.
	ckPos, ckVal, ckFound := cfFindCard(hdr, header.KeyChecksum, 0)
	if ckFound && ckVal == "" {
		return 0, 0, fmt.Errorf("%w: %s (cfitsio: VALUE_UNDEFINED)", ErrUndefinedChecksum, header.KeyChecksum)
	}
	dsFrom := 0
	if ckFound {
		dsFrom = ckPos + 1
	}
	_, dsVal, dsFound := cfFindCard(hdr, header.KeyDatasum, dsFrom)
	if dsFound && dsVal == "" {
		return 0, 0, fmt.Errorf("%w: %s (cfitsio: VALUE_UNDEFINED)", ErrUndefinedChecksum, header.KeyDatasum)
	}
	if dsFound {
		if s := cfUnquote(dsVal); s != "" {
			datasum = -1
			if f := math.Trunc(cfAtof(s)); f >= 0 && f < 0x1p32 && uint32(f) == dataSum {
				datasum = 1
			}
		}
	}
	if ckFound && cfUnquote(ckVal) != "" {
		hdusum = -1
		if s := checksum.Update(dataSum, hdr); s == 0xFFFFFFFF || s == 0 {
			hdusum = 1
		}
	}
	return datasum, hdusum, nil
}

// VerifyChecksum verifies the CHECKSUM and DATASUM keywords on an HDU if
// they are present. It returns nil if every present keyword verifies (or
// if both are absent) and an error naming the keyword(s) that do not.
// The result agrees with cfitsio's fits_verify_chksum (see ChecksumStatus).
func VerifyChecksum(h HDU) error {
	ds, hs, err := ChecksumStatus(h)
	if err != nil {
		return err
	}
	switch {
	case ds < 0 && hs < 0:
		return errors.New("fits: DATASUM and CHECKSUM do not match the HDU")
	case ds < 0:
		return errors.New("fits: DATASUM does not match the data unit")
	case hs < 0:
		return errors.New("fits: CHECKSUM does not match the HDU")
	}
	return nil
}

// onDiskSums reads the HDU's header blocks as stored and returns them with
// the ones'-complement sum of the data unit (to the end of its padding;
// a missing final pad reads as the zeros it stands for).
func (r *hduRecord) onDiskSums() (hdr []byte, dataSum uint32, err error) {
	hdr = make([]byte, r.dataStart-r.headerStart)
	if err := r.file.br.ReadRange(r.headerStart, hdr); err != nil {
		return nil, 0, fmt.Errorf("fits: checksum: read header of HDU %d: %w", r.index, err)
	}
	chunk := make([]byte, 64*block.Size)
	for off := r.dataStart; off < r.paddedEnd; {
		n := min(int64(len(chunk)), r.paddedEnd-off)
		if err := r.file.br.ReadRange(off, chunk[:n]); err != nil {
			return nil, 0, fmt.Errorf("fits: checksum: read data of HDU %d: %w", r.index, err)
		}
		dataSum = checksum.Update(dataSum, chunk[:n])
		off += n
	}
	return hdr, dataSum, nil
}

// WriteChecksum sets the DATASUM and CHECKSUM cards of an HDU of a file
// opened with OpenForEdit, the way cfitsio's fits_write_chksum computes
// them: DATASUM from the data unit on disk, CHECKSUM so that the HDU,
// with the header exactly as Flush will write it, sums to all ones. Call
// (*File).Flush afterwards to write the header; data changes made after
// WriteChecksum invalidate it. The computed DATASUM and the 32-bit value
// encoded in CHECKSUM are returned.
//
// WriteChecksum returns ErrReadOnly in ModeRead and an error in
// ModeCreate (Flush does not rewrite headers of a file being created, so
// the cards would never reach the disk): reopen the file with
// OpenForEdit to checksum it.
func WriteChecksum(h HDU) (dataSum, checkSum uint32, err error) {
	rec := hduRecordOf(h)
	if rec == nil || rec.file == nil {
		return 0, 0, errors.New("fits: WriteChecksum: HDU has no backing record")
	}
	switch rec.file.mode {
	case ModeRead:
		return 0, 0, ErrReadOnly
	case ModeCreate:
		return 0, 0, errors.New("fits: WriteChecksum: not supported in ModeCreate; close the file and reopen it with OpenForEdit")
	}
	hdr, err := rec.loadHeader()
	if err != nil {
		return 0, 0, err
	}
	_, dataSum, err = rec.onDiskSums()
	if err != nil {
		return 0, 0, err
	}
	if err := hdr.Set(header.KeyDatasum, strconv.FormatUint(uint64(dataSum), 10), "data unit checksum"); err != nil {
		return 0, 0, err
	}
	if err := hdr.Set(header.KeyChecksum, "0000000000000000", "HDU checksum"); err != nil {
		return 0, 0, err
	}
	hdrBytes, err := header.Encode(hdr)
	if err != nil {
		return 0, 0, err
	}
	// The encoding of the complement, written over the all-'0'
	// placeholder, adds exactly the complement to the sum (Appendix J).
	sum := checksum.Update(dataSum, hdrBytes)
	checkSum = ^sum
	enc := checksum.Encode(checkSum, false)
	if err := hdr.Set(header.KeyChecksum, string(enc[:]), "HDU checksum"); err != nil {
		return 0, 0, err
	}
	rec.dirty = true
	return dataSum, checkSum, nil
}

// hduRecordOf extracts the private *hduRecord from an HDU, regardless of
// concrete type. Returns nil if h is not backed by a record.
func hduRecordOf(h HDU) *hduRecord {
	switch v := h.(type) {
	case *ImageHDU:
		return v.rec
	case *BinaryTableHDU:
		return v.rec
	case *ASCIITableHDU:
		return v.rec
	case *CompressedImageHDU:
		return v.rec
	}
	return nil
}
