package main

import (
	"bufio"
    "os"
    "io"
    "bytes"
    "gopkg.in/cheggaaa/pb.v1"
)

type Read struct {
    name string
    seq string
}

type Parser struct {
    file *os.File
    scanner *bufio.Scanner
    fastq bool
    bar *pb.ProgressBar
}

func (p *Parser) GetRead() (Read, bool) {
    if !p.scanner.Scan() {
        p.bar.FinishPrint("The End!")
        return Read{"",""}, true
    }
    header := p.scanner.Text()
    p.scanner.Scan()
    sequence := p.scanner.Text()    
    if p.fastq {
        p.scanner.Scan()
        p.scanner.Scan()
    }
    p.bar.Increment()
    return Read{header[1:],sequence}, false
}

func (p *Parser) Close() {
    p.file.Close()
}

func CreateParser(referenceFile string, fastq bool) Parser {
	f, err := os.Open(referenceFile)
	if err != nil {
		panic(err)
	}
    numLines,_  := wc(f)
    f.Close()
    f, _ = os.Open(referenceFile)

	var bar *pb.ProgressBar
    if fastq {
        bar = pb.StartNew(numLines / 4)
    } else {
        bar = pb.StartNew(numLines / 2)
    }

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanLines)
    return Parser{f,scanner,fastq,bar}
}

func wc(r io.Reader) (int, error) {
	buf := make([]byte, 32*1024)
	count := 0
	lineSep := []byte{'\n'}

	for {
		c, err := r.Read(buf)
		count += bytes.Count(buf[:c], lineSep)

		switch {
		case err == io.EOF:
			return count, nil

		case err != nil:
			return count, err
		}
	}
}