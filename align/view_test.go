package align

import (
	"fmt"
	"github.com/craiglowe/gonomics/dna"
	"testing"
)

var seqA, _ = dna.StringToBases("ACGT")
var seqB, _ = dna.StringToBases("ACGT")
var seqC, _ = dna.StringToBases("ACCGT")
var seqD, _ = dna.StringToBases("ACGGT")
var seqE, _ = dna.StringToBases("AGGAGGTTTACACGTG")
var seqF, _ = dna.StringToBases("AGGAGGACACGAAAATG")
var seqG, _ = dna.StringToBases("ACTTTTTTGGT")
var seqH, _ = dna.StringToBases("AC")
var seqI, _ = dna.StringToBases("ACC")
var alignTests = []struct {
	seqOne []dna.Base
	seqTwo []dna.Base
}{
	{seqA, seqB},
	{seqA, seqC},
	{seqH, seqI},
	{seqA, seqD},
	{seqE, seqF},
	{seqD, seqG},
	{seqG, seqD},
}

func TestConstGap(t *testing.T) {
	for _, test := range alignTests {
		score, cigar := AlignConstGap(test.seqOne, test.seqTwo, defaultScores(), -430)
		prettyAlignment := View(test.seqOne, test.seqTwo, cigar)
		fmt.Printf("score=%d\n%s\n", score, prettyAlignment)
	}
}

func TestAffineGap(t *testing.T) {
	for _, test := range alignTests {
		score, cigar := AlignAffineGap(test.seqOne, test.seqTwo, defaultScores(), -400, -30)
		prettyAlignment := View(test.seqOne, test.seqTwo, cigar)
		fmt.Printf("score=%d\n%s\n", score, prettyAlignment)
	}
}
