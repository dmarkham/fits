package fits

import (
	"fmt"
	"io"
)

// CopyTo walks every HDU in f and writes its header, data and padding bytes
// to w in file order. The output is the input up to the end of the last
// HDU, so for a file with no bytes after its last HDU (TrailingBytes() ==
// 0) and a fully padded last data unit it is byte-for-byte identical to the
// input (TestCopyToByteIdentical, every committed sample; cmd/fitscopy is a
// command-line wrapper around this check). Otherwise it differs in exactly
// two ways (TestCopyToDropsTrailingBytesAndPads): bytes after the last HDU
// are not copied, and a last data unit that lacks its final padding is
// written with it, as zeros (FITS 4.0, section 3.3.2).
//
// CopyTo does not decode, re-encode or modify anything: each header is the
// bytes on disk (header changes made through Header() are included only
// after Flush), and data passes through unread. Callers that want a
// structural rebuild must use EditFile instead.
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
