package align

import (
	"fmt"
	"github.com/craiglowe/gonomics/dna"
	"github.com/craiglowe/gonomics/fasta"
	"image/color"
	"sort"
)

type keyValue struct {
	Key   string
	Value int
}

func incrementOrAdd(list []keyValue, needle string) []keyValue {
	for i, _ := range list {
		if list[i].Key == needle {
			list[i].Value++
			return list
		}
	}
	return append(list, keyValue{Key: needle, Value: 1})
}

func determineChunkColors(aln []fasta.Fasta, chunkSize int, palette color.Palette) (map[string]color.Color, error) {
	answer := make(map[string]color.Color, 0)
	list := make([]keyValue, 0)

	for i, _ := range aln {
		if len(aln[i].Seq)%chunkSize != 0 {
			return nil, fmt.Errorf("The %s sequence has a length of %d, which is not divisible by a chunkSize of %d\n", aln[i].Name, len(aln[i].Seq), chunkSize)
		}
		for chunkStart := 0; chunkStart < len(aln[i].Seq); chunkStart += chunkSize {
			chunkText := dna.BasesToString(aln[i].Seq[chunkStart:(chunkStart + chunkSize)])
			incrementOrAdd(list, chunkText)
		}
	}
	sort.Slice(list, func(i, j int) bool { return list[i].Value < list[j].Value })

	for i := 0; i < len(list) && i < len(palette); i++ {
		answer[list[i].Key] = palette[i]
	}

	return answer, nil
}
