package main

import (
	"bufio"
	"fmt"
	"os"
)

func main() {
	f, err := os.Open(os.Args[1])
	if err != nil {
		panic(err)
	}
	defer f.Close()

	reference := ""

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanLines)
	scanner.Scan()

	for scanner.Scan() {
		reference += scanner.Text()
	}
	reference += "$"

	mapper := CreateMapper(reference)

	// for _, search := range stringsToSearch {
	// 	ints := mapper.findOccurances(search)
	// 	for _, val := range ints {
	// 		fmt.Print(val, " ")
	// 	}
	// }
	// fmt.Print("\n")
}
