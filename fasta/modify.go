package fasta

import (
	"github.com/craiglowe/gonomics/dna"
)

func ReverseComplement(record Fasta) {
	record.name = fmt.Sprintf("%s_revComp", record.name)
	record.seq = dna.ReverseComplement(record.seq)
}

func ReverseComplementAll(records []Fasta) {
	
