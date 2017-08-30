package fasta

import (
	"fmt"
	"github.com/craiglowe/gonomics/dna"
)

func ReverseComplement(record *Fasta) {
	record.name = fmt.Sprintf("%s_revComp", record.name)
	dna.ReverseComplement(record.seq)
}

func ReverseComplementAll(records []Fasta) {
	for idx, _ := range records {
		ReverseComplement(&records[idx])
	}
}
