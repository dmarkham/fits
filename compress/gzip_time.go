package compress

import "time"

// zeroTime is a fixed zero time used for gzip header mtime. Using a
// constant instead of time.Now() makes encoded output bit-reproducible
// across runs, which matters for caching, testing, and determinism.
var zeroTime = time.Time{}
