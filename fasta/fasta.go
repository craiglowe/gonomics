package fasta

import (
	"bufio"
	"fmt"
	"github.com/craiglowe/gonomics/common"
	"github.com/craiglowe/gonomics/dna"
	"log"
	"os"
	"strings"
)

type Fasta struct {
	Name string
	Seq  []dna.Base
}

func Read(filename string) ([]Fasta, error) {
	var line string
	var currSeq []dna.Base
	var answer []Fasta
	var seqIdx int64 = -1

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		switch {
		case strings.HasPrefix(line, "#"):
			// comment line in fasta file
		case strings.HasPrefix(line, ">"):
			answer = append(answer, Fasta{Name: line[1:len(line)]})
			seqIdx++
		default:
			currSeq, err = dna.StringToBases(line)
			if err != nil {
				log.Fatal(err)
			}
			answer[seqIdx].Seq = append(answer[seqIdx].Seq, currSeq...)
		}
	}
	return answer, scanner.Err()
}

func Write(filename string, records []Fasta) {
	file, err := os.Create(filename)
	common.ExitIfError(err)
	defer file.Close()

	for _, rec := range records {
		fmt.Fprintf(file, ">%s\n%s\n", rec.Name, dna.BasesToString(rec.Seq))
	}
}
