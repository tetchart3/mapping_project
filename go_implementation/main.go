package main

import (
	"fmt"
	"os"
	"flag"
)

func main() {

	referenceFile := flag.String("reference", "", "Reference File")
	sequenceFile := flag.String("sequence", "", "Sequence")
	fastq := flag.Bool("fastq", false, "Are the reads in Fastq format? Otherwise assumes Fasta.")
	kPointer := flag.Int("k",5,"The K to break kmers into")

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


	mapper := CreateMapper(*referenceFile)
	parser := CreateParser(*sequenceFile,*fastq)

	defer parser.Close()
	read, done := parser.GetRead();
	for !done {
		starts := make(map[int]int)
		for i := 0; i <= len(read.seq) - k; i++ {
			kmer := read.seq[i:i+k]
			indices := mapper.findOccurances(kmer)
			for _, index := range indices {
				_, exists := starts[index - i]
				if !exists {
					starts[index - i] = 0
				}
				starts[index - i]++ 
			}
		}
		fmt.Println(read.seq,starts)
		// fmt.Println(mapper.findOccurances(read.seq))
		// fmt.Println(indices)
		read, done = parser.GetRead();
	}
}
