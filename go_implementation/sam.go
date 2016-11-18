package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
)

const (
	MAPQ  = 255
	RNEXT = "*"
	PNEXT = 0
	TLEN  = 0
	QUAL  = "*"
)

type SamWriter struct {
	file    *os.File
	writer  *bufio.Writer
	channel *chan string
}

func CreateSamWriter(samFile string) SamWriter {
	f, err := os.Create(samFile)
	if err != nil {
		panic(err)
	}

	writer := bufio.NewWriter(f)
	channel := make(chan string, 500)
	return SamWriter{f, writer, &channel}
}

func (sw *SamWriter) Done() {
	close(*sw.channel)
}

func (sw *SamWriter) WriteLine(qname string, rname string, start int, qseq string, rseq string) {
	var cigar string
	var flag int
	if len(rseq) == 0 {
		cigar = "*"
		flag = 4
	} else {
		cigar = makeCigar(rseq, qseq)
		flag = 0
	}
	msg := fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
		qname, flag, rname, start+1, MAPQ, cigar, RNEXT, PNEXT, TLEN, qseq, QUAL)
	*sw.channel <- msg
}

func (sw *SamWriter) WriteSam() {
	for msg := range *sw.channel {

		_, err := sw.writer.WriteString(msg)
		if err != nil {
			panic(err)
		}

		sw.writer.Flush()

		if err != nil {
			panic(err)
		}

	}
	sw.writer.Flush()
	sw.file.Close()
}

func makeCigar(reference string, query string) string {
	cigar := ""
	currentMatchType := ""
	currentLength := 0
	for i := 0; i < len(reference); i++ {
		var matchType string
		if reference[i] == query[i] {
			matchType = "M"
		} else {
			matchType = "X"
		}
		if currentMatchType == "" {
			currentMatchType = matchType
			currentLength++
		} else if matchType == currentMatchType {
			currentLength++
		} else {
			cigar += (strconv.Itoa(currentLength) + currentMatchType)
			currentMatchType = matchType
			currentLength = 1
		}
	}
	cigar += (strconv.Itoa(currentLength) + currentMatchType)
	return cigar
}
