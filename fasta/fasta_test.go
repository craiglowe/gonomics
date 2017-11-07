package fasta

import (
	"github.com/craiglowe/gonomics/dna"
	"testing"
)

var seqOneA, _ = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqOneB, _ = dna.StringToBases("acgtACGTACGT")
var seqOneC, _ = dna.StringToBases("ACGTACGTACGTT")
var readWriteTests = []struct {
	filename string // input
	data     []Fasta
}{
	{"testdata/testOne.fa", []Fasta{{"apple", seqOneA}, {"banana", seqOneB}, {"carrot", seqOneC}}},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual, err := Read(test.filename)
		if err != nil {
			t.Errorf("Reading %s gave an error..", test.filename)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	for _, test := range readWriteTests {
		Write(test.filename+".tmp", test.data)
		actual, err := Read(test.filename + ".tmp")
		if err != nil {
			t.Errorf("Reading %s gave an error.", test.filename)
		}
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}
