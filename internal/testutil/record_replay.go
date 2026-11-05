//go:build !interop

package testutil

import "testing"

// Tier 1: a tool reading is always the committed recording.
func execute(t testing.TB, tool Tool, c *Call, r *record, _ placeholders) *record {
	t.Helper()
	return replay(t, tool, c, r)
}
