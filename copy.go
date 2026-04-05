package fits

import (
	"fmt"
	"io"
)

// CopyTo walks every HDU in f and writes its raw header + data + pad bytes
// to w in file order. When f was opened from a valid FITS file, the output
// is byte-for-byte identical to the input — this is the heart of the
// regression tool cmd/fitscopy.
//
// CopyTo does not re-encode or modify anything: headers pass through from
// rec.rawHeader (which still holds the on-disk bytes unless the caller has
// mutated the parsed header and called Flush), data passes through via
// ReadRange on the block reader. Callers that want a structural rebuild
// must use EditFile instead.
func (f *File) CopyTo(w io.Writer) error {
	for _, rec := range f.hdus {
		if _, err := w.Write(rec.rawHeader); err != nil {
			return fmt.Errorf("fits: CopyTo HDU %d header: %w", rec.index, err)
		}
		dataLen := rec.paddedEnd - rec.dataStart
		if dataLen == 0 {
			continue
		}
		const chunk = 1 << 20
		buf := make([]byte, chunk)
		var done int64
		for done < dataLen {
			n := int64(chunk)
			if dataLen-done < n {
				n = dataLen - done
			}
			if err := f.br.ReadRange(rec.dataStart+done, buf[:n]); err != nil {
				return fmt.Errorf("fits: CopyTo HDU %d data: %w", rec.index, err)
			}
			if _, err := w.Write(buf[:n]); err != nil {
				return fmt.Errorf("fits: CopyTo HDU %d write: %w", rec.index, err)
			}
			done += n
		}
	}
	return nil
}
