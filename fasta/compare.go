package fasta

import (
	"github.com/craiglowe/gonomics/dna"
	"sort"
	"strings"
)

func CompareName(alpha Fasta, beta Fasta) int {
	return strings.Compare(alpha.name, beta.name)
}

func CompareSeq(alpha Fasta, beta Fasta) int {
	return dna.CompareSeqs(alpha.seq, beta.seq, true)
}

func CompareSeqIgnoreCase(alpha Fasta, beta Fasta) int {
	return dna.CompareSeqs(alpha.seq, beta.seq, false)
}

func SortByName(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return CompareName(seqs[i], seqs[j]) == -1 })
}

func SortBySeq(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return CompareSeqIgnoreCase(seqs[i], seqs[j]) == -1 })
}
