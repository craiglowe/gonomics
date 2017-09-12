package align

import (
	"github.com/craiglowe/gonomics/common"
	"github.com/craiglowe/gonomics/dna"
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
	runLength uint64
	op        colType
}

type alignParamLinear struct {
	score [5][5]int64
	gap   int64
}

type alignParamAffine struct {
	score  [5][5]int64
	open   int64
	extend int64
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

func AlignConstGap(alpha []dna.Base, beta []dna.Base, scores [5][5]int64, gapPen int64) (int64, []cigar) {
	m := make([][]int64, len(alpha)+1)
	trace := make([][]colType, len(alpha)+1)
	for idx := range m {
		m[idx] = make([]int64, len(beta)+1)
		trace[idx] = make([]colType, len(beta)+1)
	}

	var i, j, routeIdx int
	for i, _ = range m {
		for j, _ = range m[0] {
			if i == 0 && j == 0 {
				m[i][j] = 0
			} else if i == 0 {
				m[i][j] = m[i][j-1] + gapPen
				trace[i][j] = 1
			} else if j == 0 {
				m[i][j] = m[i-1][j] + gapPen
				trace[i][j] = 2
			} else {
				m[i][j], trace[i][j] = tripleMaxTrace(m[i-1][j-1]+scores[alpha[i-1]][beta[j-1]], m[i][j-1]+gapPen, m[i-1][j]+gapPen)
			}
		}
	}

	route := make([]cigar, 1)
	for i, j, routeIdx = len(trace)-1, len(trace[0])-1, 0; i > 0 || j > 0; {
		if route[routeIdx].runLength == 0 {
			route[routeIdx].runLength = 1
			route[routeIdx].op = trace[i][j]
		} else if route[routeIdx].op == trace[i][j] {
			route[routeIdx].runLength += 1
		} else {
			route = append(route, cigar{runLength: 1, op: trace[i][j]})
			routeIdx++
		}
		switch trace[i][j] {
		case 0:
			i, j = i-1, j-1
		case 1:
			j -= 1
		case 2:
			i -= 1
		default:
			common.Exit("Error: unexpected traceback")
		}
	}
	reverseCigar(route)
	return m[len(m)-1][len(m[0])-1], route
}

func AlignAffineGap(alpha []dna.Base, beta []dna.Base, scores [5][5]int64, gapOpen int64, gapExtend int64) (int64, []cigar) {
	// the data structure is a 3d slice where the first index is 0,1,2 and represents
	// the match, gap in x (first seq), and gap in y (second seq).
	m := make([][][]int64, 3)
	trace := make([][][]colType, 3)
	for k, _ := range m {
		m[k] = make([][]int64, len(alpha)+1)
		trace[k] = make([][]colType, len(alpha)+1)
		for i, _ := range m[0] {
			m[k][i] = make([]int64, len(beta)+1)
			trace[k][i] = make([]colType, len(beta)+1)
		}
	}

	/*m := make([][]int64, len(alpha)+1)
		x := make([][]int64, len(alpha)+1)
		y := make([][]int64, len(alpha)+1)
	        traceM := make([][]colType, len(alpha)+1)
		traceX := make([][]colType, len(alpha)+1)
		traceY := make([][]colType, len(alpha)+1)
	        for idx := range m {
	                m[idx] = make([]int64, len(beta)+1)
			x[idx] = make([]int64, len(beta)+1)
			y[idx] = make([]int64, len(beta)+1)
	                traceM[idx] = make([]colType, len(beta)+1)
			traceX[idx] = make([]colType, len(beta)+1)
			traceY[idx] = make([]colType, len(beta)+1)
	        }*/

	for i, _ := range m[0] {
		for j, _ := range m[0][0] {
			if i == 0 && j == 0 {
				m[0][i][j] = 0
				m[1][i][j] = gapOpen
				m[2][i][j] = gapOpen
			} else if i == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = gapExtend + m[1][i][j-1]
				m[2][i][j] = veryNegNum
				trace[2][i][j] = colI
			} else if j == 0 {
				m[0][i][j] = veryNegNum
				m[1][i][j] = veryNegNum
				trace[1][i][j] = colD
				m[2][i][j] = gapExtend + m[2][i-1][j]
			} else {
				m[0][i][j], trace[0][i][j] = tripleMaxTrace(scores[alpha[i-1]][beta[j-1]]+m[0][i-1][j-1], scores[alpha[i-1]][beta[j-1]]+m[1][i-1][j-1], scores[alpha[i-1]][beta[j-1]]+m[2][i-1][j-1])
				m[1][i][j], trace[1][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i][j-1], gapExtend+m[1][i][j-1], gapOpen+gapExtend+m[2][i][j-1])
				m[2][i][j], trace[2][i][j] = tripleMaxTrace(gapOpen+gapExtend+m[0][i-1][j], gapOpen+gapExtend+m[1][i-1][j], gapExtend+m[2][i-1][j])
			}
		}
	}

	route := make([]cigar, 1)
	lastI := len(m[0]) - 1
	lastJ := len(m[0][0]) - 1
	maxScore, k := tripleMaxTrace(m[0][lastI][lastJ], m[1][lastI][lastJ], m[2][lastI][lastJ])
	for i, j, routeIdx := lastI, lastJ, 0; i > 0 || j > 0; {
		if route[routeIdx].runLength == 0 {
			route[routeIdx].runLength = 1
			route[routeIdx].op = k
		} else if route[routeIdx].op == k {
			route[routeIdx].runLength += 1
		} else {
			route = append(route, cigar{runLength: 1, op: k})
			routeIdx++
		}
		switch k {
		case colM:
			k = trace[k][i][j]
			i--
			j--
		case colI:
			k = trace[k][i][j]
			j--
		case colD:
			k = trace[k][i][j]
			i--
		default:
			common.Exit("Error: unexpected traceback")
		}
	}

	reverseCigar(route)
	return maxScore, route
}
