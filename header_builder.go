package fits

import "github.com/dmarkham/fits/header"

// headerBuilder builds an HDU header and keeps the first error of its
// Set / AddCard calls, so a writer can make a run of calls and check once.
// Every writer must check err before encoding: a card that was refused is
// otherwise missing from the HDU without a word.
type headerBuilder struct {
	h   *header.Header
	err error
}

func newHeaderBuilder() *headerBuilder { return &headerBuilder{h: header.New()} }

func (b *headerBuilder) set(name string, value any, comment string) {
	if b.err == nil {
		b.err = b.h.Set(name, value, comment)
	}
}

// add appends a caller-supplied card. A card read from a file keeps its
// original bytes (header.Header.AddCard); one built in memory is validated.
//
// A HIERARCH alias of a standard keyword is refused with the rule Flush
// applies (refuseHierarchAlias), so no writer emits a card that cfitsio
// and astropy read differently.
func (b *headerBuilder) add(c header.Card) {
	if b.err == nil {
		b.err = refuseHierarchAlias(c.Key)
	}
	if b.err == nil {
		b.err = b.h.AddCard(c)
	}
}
