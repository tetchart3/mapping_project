package main

import (
	"flag"
	"fmt"
	"os"
	"runtime"
	"sync"
)

const (
	AGREEMENT = 10
)

func main() {

	referenceFile := flag.String("reference", "", "Reference File")
	sequenceFile := flag.String("reads", "", "Sequence")
	samFile := flag.String("sam", "", "SAM File to write")
	fastq := flag.Bool("fastq", false, "Are the reads in Fastq format? Otherwise assumes Fasta.")
	kPointer := flag.Int("k", 5, "The K to break kmers into")

	flag.Parse()

	var k int
	k = *kPointer

	if len(*referenceFile) == 0 {
		fmt.Println("Reference file required.")
		os.Exit(1)
	}

	if len(*sequenceFile) == 0 {
		fmt.Println("Sequence file required.")
		os.Exit(1)
	}

	if len(*samFile) == 0 {
		fmt.Println("SAM file required.")
		os.Exit(1)
	}

	mapper := CreateMapper(*referenceFile)
	parser := CreateParser(*sequenceFile, *fastq)
	writer := CreateSamWriter(*samFile)
	var wg sync.WaitGroup

	go writer.WriteSam()
	unmappedReads := 0
	reverseReads := 0
	avgAgreement := 0
	i := 0
	read, done := parser.GetRead()
	for !done {
		wg.Add(1)
		currentRead := Read{read.name, read.seq}
		go analyzeRead(&mapper, &writer, &currentRead, k, &unmappedReads, &reverseReads, &avgAgreement, &wg)

		if (i % (runtime.NumCPU() * 2)) == 0 {
			wg.Wait()
		}
		i++
		read, done = parser.GetRead()
	}
	wg.Wait()
	parser.Close()
	writer.Done()
	fmt.Println(unmappedReads, "unmapped reads.")
	fmt.Println(reverseReads, "reverse reads.")
	fmt.Println(avgAgreement/i, "average agreement")
}

func analyzeRead(mapper *Mapper, writer *SamWriter, read *Read, k int, unmappedReads *int, reverseReads *int, avgAgreement *int, wg *sync.WaitGroup) {
	defer wg.Done()
	var bestStart int
	forwardStart, forwardAgreement := findStart(mapper, *read, k)
	if forwardAgreement < (len(read.seq)-k)/3 {
		revRead := read.ReverseCompliment()
		reverseStart, reverseAgreement := findStart(mapper, revRead, k)

		if forwardAgreement > reverseAgreement {
			bestStart = forwardStart
			*avgAgreement += forwardAgreement
		} else {
			bestStart = reverseStart
			*reverseReads++
			*avgAgreement += reverseAgreement

		}
	} else {
		bestStart = forwardStart
		*avgAgreement += forwardAgreement
	}

	if bestStart >= 0 && (bestStart+len(read.seq)) < len(mapper.reference)-1 {
		writer.WriteLine(read.name, mapper.referenceName, bestStart, read.seq, mapper.reference[bestStart:bestStart+len(read.seq)])
	} else {
		writer.WriteLine(read.name, "*", -1, read.seq, "")
		*unmappedReads++
	}
}

func findStart(mapper *Mapper, read Read, k int) (int, int) {
	starts := make(map[int]int)
	for i := 0; i <= len(read.seq)-k; i++ {
		kmer := read.seq[i : i+k]
		indices := mapper.findOccurances(kmer)
		for _, index := range indices {
			_, exists := starts[index-i]
			if !exists {
				starts[index-i] = 0
			}
			starts[index-i]++
		}
	}
	// fmt.Println(starts)
	bestStart := -1
	bestAgreement := 0
	for start, agreement := range starts {
		if agreement > bestAgreement {
			bestStart = start
			bestAgreement = agreement
		}
	}
	// if bestAgreement != 46 {
	// 	fmt.Println(len(read.seq), bestStart, bestAgreement)
	// }
	return bestStart, bestAgreement
}
