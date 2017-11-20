package align

import (
	"math"
)

var veryNegNum int64 = math.MinInt64 / 2

// these are relative to the first seq.
// e.g. colI is an insertion in the second seq, relative to the first
type colType uint8

const (
	colM colType = 0
	colI colType = 1
	colD colType = 2
)

type cigar struct {
	runLength int64
	op        colType
}

func defaultScores() [5][5]int64 {
	answer := [5][5]int64{
		{91, -114, -31, -123, -44},
		{-114, 100, -125, -31, -43},
		{-31, -125, 100, -114, -43},
		{-123, -31, -114, 91, -44},
		{-44, -43, -43, -44, -43},
	}
	return answer
}

func tripleMaxTrace(a int64, b int64, c int64) (int64, colType) {
	if a >= b && a >= c {
		return a, colM
	} else if b >= c {
		return b, colI
	} else {
		return c, colD
	}
}

func reverseCigar(alpha []cigar) {
	for i, j := 0, len(alpha)-1; i < j; i, j = i+1, j-1 {
		alpha[i], alpha[j] = alpha[j], alpha[i]
	}
}

func countAlignmentColumns(route []cigar) int64 {
	var count int64 = 0
	for i, _ := range route {
		count += route[i].runLength
	}
	return count
}
