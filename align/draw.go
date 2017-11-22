package align

import (
	"fmt"
	"github.com/craiglowe/gonomics/dna"
	"github.com/craiglowe/gonomics/draw"
	"github.com/craiglowe/gonomics/fasta"
	"image"
	"image/color"
	"sort"
	"strings"
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
			gapCount := strings.Count(chunkText, "-")
			switch gapCount {
				case chunkSize: /* ignore if all gaps */
				case 0:
					list = incrementOrAdd(list, chunkText)
				default:
					return nil, fmt.Errorf("Error: %s should be either all gaps or no gaps\n", chunkText)
			}
		}
	}
	fmt.Printf("Number of chunks recorded: %d\n", len(list))
	sort.Slice(list, func(i, j int) bool { return list[i].Value > list[j].Value })

	for i := 0; i < len(list) && i < len(palette); i++ {
		answer[list[i].Key] = palette[i]
	}

	return answer, nil
}

func DrawAlignedChunks(aln []fasta.Fasta, chunkSize int, chunkPixelWidth int, chunkPixelHeight int, border int, seqPixelSpacing int) (*image.RGBA, error) {
	colorMap, err := determineChunkColors(aln, chunkSize, draw.TrubetskoyPalette[:19])
	fmt.Printf("colors to be used: %v\n", colorMap)
	if err != nil {
		return nil, err
	}
	allGaps := strings.Repeat("-", chunkSize)
	colorMap[allGaps] = draw.TrubetskoyPalette[21] /* black */

	alnLength := len(aln[0].Seq)
	numSeq := len(aln)
	imageWidth := border*2 + alnLength / chunkSize * chunkPixelWidth
	imageHeight := border*2 + chunkPixelHeight * numSeq + seqPixelSpacing * (numSeq-1)
	img := image.NewRGBA(image.Rect(0, 0, imageWidth, imageHeight))
	draw.FilledRectangle(img, 0, 0, imageWidth, imageHeight, draw.TrubetskoyPalette[20]) /* make everything white to start */
	for i, _ := range aln {
		for chunkStart := 0; chunkStart < len(aln[i].Seq); chunkStart += chunkSize {
			chunkText := dna.BasesToString(aln[i].Seq[chunkStart:(chunkStart + chunkSize)])
			chunkColor, found := colorMap[chunkText]
			if !found {
				chunkColor = draw.TrubetskoyPalette[19] /* gray */
			}
			xStart := border + chunkStart/chunkSize*chunkPixelWidth
			xEnd := xStart + chunkPixelWidth
			yStart := border + i*chunkPixelHeight + (i-1)*seqPixelSpacing
			yEnd := yStart + chunkPixelHeight
			draw.FilledRectangle(img, xStart, yStart, xEnd, yEnd, chunkColor)
		}
	}
	return img, nil
}
