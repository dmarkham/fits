//go:build !linux && !darwin

package fits

// copyXattrs is not implemented outside Linux and macOS, where editing is
// refused (editSupported); it refuses too rather than drop attributes
// silently.
func copyXattrs(src, dst string) error { return editSupported() }
