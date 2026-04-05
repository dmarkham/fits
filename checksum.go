package fits

import (
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/checksum"
)

// VerifyChecksum verifies the CHECKSUM and DATASUM keywords on an HDU if
// they are present. It returns nil if both verify (or if they are absent).
// If one or both are present but wrong, it returns an error identifying
// which failed.
//
// Algorithm (FITS Appendix J):
//   - DATASUM is the 32-bit 1's-complement checksum of the data section
//     only, encoded as an ASCII integer.
//   - CHECKSUM is the 16-char ASCII-encoded 1's-complement checksum such
//     that the combined (header + data) checksum sums to all-ones.
func VerifyChecksum(h HDU) error {
	rec := hduRecordOf(h)
	if rec == nil {
		return errors.New("fits: VerifyChecksum: HDU has no backing record")
	}
	hdr := h.Header()
	var dataSumVal uint32
	var haveDataSum bool
	if s, err := hdr.String(header.KeyDatasum); err == nil {
		haveDataSum = true
		s = strings.TrimSpace(s)
		// cfitsio writes DATASUM as a decimal integer string.
		v, convErr := strconv.ParseUint(s, 10, 32)
		if convErr != nil {
			return fmt.Errorf("fits: DATASUM not a valid integer: %q", s)
		}
		dataSumVal = uint32(v)
	}

	// Read the data section in 2880-byte chunks and accumulate.
	var actualDataSum uint32
	dataBytes := rec.paddedEnd - rec.dataStart
	if dataBytes > 0 {
		chunk := make([]byte, 2880)
		off := rec.dataStart
		for off < rec.paddedEnd {
			if err := rec.file.br.ReadRange(off, chunk); err != nil {
				return fmt.Errorf("fits: VerifyChecksum: read data: %w", err)
			}
			actualDataSum = checksum.Update(actualDataSum, chunk)
			off += 2880
		}
	}

	if haveDataSum && actualDataSum != dataSumVal {
		return fmt.Errorf("fits: DATASUM mismatch: stored %d, computed %d", dataSumVal, actualDataSum)
	}

	// CHECKSUM is the sum over header bytes + data bytes such that the
	// result is complementary. If present, verify by summing header bytes
	// and checking that header_sum + data_sum + stored-card = all-ones.
	checkStr, err := hdr.String(header.KeyChecksum)
	if err != nil {
		return nil // no CHECKSUM card; DATASUM already verified (or also absent)
	}
	if len(checkStr) != 16 {
		return fmt.Errorf("fits: CHECKSUM length %d != 16", len(checkStr))
	}
	// Sum the header bytes (as currently on disk).
	hdrBytes := rec.rawHeader
	headerSum := checksum.Compute(hdrBytes)
	// The CHECKSUM card itself was part of the header during the compute —
	// cfitsio uses a "placeholder trick": write "0000...0000" first,
	// compute, then derive the final CHECKSUM value so that the sum over
	// the entire HDU is zero (complement). On verify we read back the
	// stored CHECKSUM, decode its value, and check that
	// header_sum + data_sum + stored_sum == 0xFFFFFFFF.
	var enc [16]byte
	copy(enc[:], checkStr)
	storedSum := checksum.Decode(enc, false)
	total := headerSum + actualDataSum + storedSum
	if total != 0xFFFFFFFF && total != 0 {
		return fmt.Errorf("fits: CHECKSUM mismatch: combined sum %#x", total)
	}
	return nil
}

// WriteChecksum updates (or creates) the CHECKSUM and DATASUM cards on an
// HDU in edit mode. It must be called BEFORE Flush so that the re-serialized
// header contains the correct cards. The computed values are returned for
// inspection.
//
// Checksums are not computed automatically on Close — callers must opt
// in via WriteChecksum per HDU. This keeps the strict-fidelity promise:
// users get byte-exact round-trip unless they explicitly ask for mutation.
func WriteChecksum(h HDU) (dataSum, checkSum uint32, err error) {
	rec := hduRecordOf(h)
	if rec == nil {
		return 0, 0, errors.New("fits: WriteChecksum: HDU has no backing record")
	}
	if rec.file == nil || rec.file.mode == ModeRead {
		return 0, 0, ErrReadOnly
	}
	hdr := h.Header()

	// Compute DATASUM by reading the on-disk data bytes.
	dataBytes := rec.paddedEnd - rec.dataStart
	if dataBytes > 0 {
		chunk := make([]byte, 2880)
		off := rec.dataStart
		for off < rec.paddedEnd {
			if err := rec.file.br.ReadRange(off, chunk); err != nil {
				return 0, 0, err
			}
			dataSum = checksum.Update(dataSum, chunk)
			off += 2880
		}
	}

	// Update DATASUM card.
	dsStr := strconv.FormatUint(uint64(dataSum), 10)
	if err := hdr.Set(header.KeyDatasum, dsStr, "data checksum"); err != nil {
		return 0, 0, err
	}

	// Placeholder CHECKSUM and re-encode the header to compute header sum.
	if err := hdr.Set(header.KeyChecksum, "0000000000000000", "HDU checksum"); err != nil {
		return 0, 0, err
	}
	hdrBytes, err := header.Encode(hdr)
	if err != nil {
		return 0, 0, err
	}
	headerSum := checksum.Compute(hdrBytes)
	total := headerSum + dataSum
	// We want total + stored == 0xFFFFFFFF; stored = 0xFFFFFFFF - total.
	stored := 0xFFFFFFFF - total
	checkSum = stored
	checkEnc := checksum.Encode(stored, false)
	if err := hdr.Set(header.KeyChecksum, string(checkEnc[:]), "HDU checksum"); err != nil {
		return 0, 0, err
	}
	// Mark dirty so a subsequent Flush persists the updated header.
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
	}
	return nil
}
