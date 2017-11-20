package main

import (
	"flag"
	"fmt"
	"github.com/craiglowe/gonomics/fasta"
	"github.com/craiglowe/gonomics/align"
	"log"
	"strconv"
	"github.com/craiglowe/gonomics/common"
)

func faChunkAlign(inFile string, chunkSize int, outFile string) {
	records, err := fasta.Read(inFile)
	common.ExitIfError(err)
	fmt.Printf("len = %d\n", len(records))

	records = align.AllSeqAffineChunk(records, chunkSize)
	fmt.Printf("len = %d\n", len(records))

	fasta.Write(outFile, records)
}

func usage() {
	fmt.Print(
		"faChunkAlign - Align two or more sequeces by \"chunks\" of bases\n" +
		"                instead of by single bases.  Each sequence must\n" +
		"                have a length that is divisible by the chunk size.\n" +
		"Usage:\n" +
		" faChunkAlign unaligned.fa chunkSize aligned.fa\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chunkSize, err := strconv.Atoi(flag.Arg(1))
	common.ExitIfError(err)
	outFile := flag.Arg(2)

	faChunkAlign(inFile, chunkSize, outFile)
}
