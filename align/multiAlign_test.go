package align

import (
	"github.com/craiglowe/gonomics/fasta"
	"testing"
)

var multiAlignTests = []struct {
	input    string
	expected string
}{
	{"testdata/multiAlignTest.in.fa", "testdata/multiAlignTest.expected.fa"},
}

func TestMultiAlignGap(t *testing.T) {
	for _, test := range multiAlignTests {
		input, err := fasta.Read(test.input)
		if err != nil {
			t.Errorf("Reading %s gave an error..", test.input)
		}
		expected, err := fasta.Read(test.expected)
		if err != nil {
			t.Errorf("Reading %s gave an error..", test.expected)
		}

		aligned := AllSeqAffine(input)

		if !fasta.AllAreEqualIgnoreOrder(aligned, expected) {
			fasta.Write("testdata/multiAlignTest.tmp", aligned)
			t.Error("Alignment not as expected: testdata/multiAlignTest.tmp does not equal %s.", test.expected)
		}
	}
}
