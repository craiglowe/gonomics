package fasta

import (
	"fmt"
	"github.com/craiglowe/gonomics/dna"
)

func ReverseComplement(record *Fasta) {
	record.Name = fmt.Sprintf("%s_revComp", record.Name)
	dna.ReverseComplement(record.Seq)
}

func ReverseComplementAll(records []Fasta) {
	for idx, _ := range records {
		ReverseComplement(&records[idx])
	}
}
