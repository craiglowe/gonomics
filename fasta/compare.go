package fasta

import (
	"github.com/craiglowe/gonomics/dna"
	"sort"
	"strings"
)

func compareName(alpha Fasta, beta Fasta) int {
	return strings.Compare(alpha.name, beta.name)
}

func compareSeq(alpha Fasta, beta Fasta) int {
	return dna.CompareSeqsCaseSensitive(alpha.seq, beta.seq)
}

func compareSeqIgnoreCase(alpha Fasta, beta Fasta) int {
	return dna.CompareSeqsIgnoreCase(alpha.seq, beta.seq)
}

func isEqual(alpha Fasta, beta Fasta) bool {
	if compareName(alpha, beta) == 0 && compareSeq(alpha, beta) == 0 {
		return true
	} else {
		return false
	}
}

func AllAreEqual(alpha []Fasta, beta []Fasta) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for idx, _ := range alpha {
		if !isEqual(alpha[idx], beta[idx]) {
			return false
		}
	}
	return true
}

func SortByName(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareName(seqs[i], seqs[j]) == -1 })
}

func SortBySeq(seqs []Fasta) {
	sort.Slice(seqs, func(i, j int) bool { return compareSeqIgnoreCase(seqs[i], seqs[j]) == -1 })
}
