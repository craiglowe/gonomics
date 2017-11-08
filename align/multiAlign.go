package align

import (
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
	bestScore = 0
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

func AllSeqAffine(records []fasta.Fasta) []fasta.Fasta {
	groups := fastaListToIndividualGroups(records)
	for len(groups) > 1 {
		x, y, _, route := nearestGroups(groups)
		groups = mergeFastaGroups(groups, x, y, route)
	}
	return groups[0]
}
