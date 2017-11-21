package align

import (
	"math"
	"github.com/craiglowe/gonomics/dna"
	"github.com/craiglowe/gonomics/fasta"
)

func fastaListToIndividualGroups(records []fasta.Fasta) [][]fasta.Fasta {
	answer := make([][]fasta.Fasta, len(records))
	for i, _ := range answer {
		answer[i] = make([]fasta.Fasta, 1)
		answer[i][0] = records[i]
	}
	return answer
}

func mergeFastaGroups(groups [][]fasta.Fasta, x int, y int, route []cigar) [][]fasta.Fasta {
	groups[x] = mergeMultipleAlignments(groups[x], groups[y], route)
	groups[y] = groups[len(groups)-1]
	groups = groups[:len(groups)-1]
	return groups
}

func nearestGroups(groups [][]fasta.Fasta) (bestX int, bestY int, bestScore int64, bestRoute []cigar) {
	var route []cigar
	var score int64
	bestScore = math.MinInt64
	for x := 0; x < len(groups)-1; x++ {
		for y := x + 1; y < len(groups); y++ {
			score, route = multipleAffineGap(groups[x], groups[y], defaultScores(), -400, -30)
			if score > bestScore {
				bestX, bestY, bestScore, bestRoute = x, y, score, route
			}
		}
	}
	return bestX, bestY, bestScore, bestRoute
}

func nearestGroupsChunk(groups [][]fasta.Fasta, chunkSize int) (bestX int, bestY int, bestScore int64, bestRoute []cigar) {
	var route []cigar
	var score int64
	bestScore = math.MinInt64
	for x := 0; x < len(groups)-1; x++ {
		for y := x + 1; y < len(groups); y++ {
			score, route = multipleAffineGapChunk(groups[x], groups[y], defaultScores(), -400, -30, int64(chunkSize))
			if score > bestScore {
				bestX, bestY, bestScore, bestRoute = x, y, score, route
			}
		}
	}
	return bestX, bestY, bestScore, bestRoute
}

func AllSeqAffine(records []fasta.Fasta) []fasta.Fasta {
	groups := fastaListToIndividualGroups(records)
	for len(groups) > 1 {
		x, y, _, route := nearestGroups(groups)
		groups = mergeFastaGroups(groups, x, y, route)
	}
	return groups[0]
}

func AllSeqAffineChunk(records []fasta.Fasta, chunkSize int) []fasta.Fasta {
	groups := fastaListToIndividualGroups(records)
	for len(groups) > 1 {
		x, y, _, route := nearestGroupsChunk(groups, chunkSize)
		groups = mergeFastaGroups(groups, x, y, route)
	}
	return groups[0]
}

//average of pairs scoring scheme where gaps are ignored
//maybe there should be a small penalty for gaps so that gaps will tend to be in the same location
func scoreColumnMatch(alpha []fasta.Fasta, beta []fasta.Fasta, alphaCol int, betaCol int, scores [5][5]int64) int64 {
	var sum, count int64 = 0, 0
	for alphaSeqIdx, _ := range alpha {
		for betaSeqIdx, _ := range beta {
			if alpha[alphaSeqIdx].Seq[alphaCol] != dna.Gap && beta[betaSeqIdx].Seq[betaCol] != dna.Gap {
				sum += scores[alpha[alphaSeqIdx].Seq[alphaCol]][beta[betaSeqIdx].Seq[betaCol]]
				count++
			}
		}
	}
	return sum / count
}

func ungappedRegionColumnScore(alpha []fasta.Fasta, alphaStart int, beta []fasta.Fasta, betaStart int, length int, scores [5][5]int64) int64 {
	var answer int64 = 0
	for i, j := alphaStart, betaStart; i < alphaStart+length; i, j = i+1, j+1 {
		answer += scoreColumnMatch(alpha, beta, i, j, scores)
	}
	return answer
}

func mergeMultipleAlignments(alpha []fasta.Fasta, beta []fasta.Fasta, route []cigar) []fasta.Fasta {
	answer := make([]fasta.Fasta, len(alpha)+len(beta))
	totalCols := countAlignmentColumns(route)

	for i, _ := range answer {
		if i < len(alpha) {
			answer[i] = fasta.Fasta{Name: alpha[i].Name, Seq: make([]dna.Base, totalCols)}
		} else {
			answer[i] = fasta.Fasta{Name: beta[i-len(alpha)].Name, Seq: make([]dna.Base, totalCols)}
		}
	}

	var alphaCol, betaCol, ansCol int = 0, 0, 0
	for i, _ := range route {
		for j := 0; j < int(route[i].runLength); j++ {
			for k, _ := range answer {
				if k < len(alpha) {
					if route[i].op == colM || route[i].op == colD {
						answer[k].Seq[ansCol] = alpha[k].Seq[alphaCol]
					} else {
						answer[k].Seq[ansCol] = dna.Gap
					}
				} else {
					if route[i].op == colM || route[i].op == colI {
						answer[k].Seq[ansCol] = beta[k-len(alpha)].Seq[betaCol]
					} else {
						answer[k].Seq[ansCol] = dna.Gap
					}
				}
			}
			switch route[i].op {
			case colM:
				alphaCol, betaCol, ansCol = alphaCol+1, betaCol+1, ansCol+1
			case colI:
				betaCol, ansCol = betaCol+1, ansCol+1
			case colD:
				alphaCol, ansCol = alphaCol+1, ansCol+1
			}
		}
	}
	return answer
}
