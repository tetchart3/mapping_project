package main

import (
	"bufio"
	"os"
	"sort"
	"strings"
)

// Copyright 2011 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// This algorithm is based on "Faster Suffix Sorting"
//   by N. Jesper Larsson and Kunihiko Sadakane
// paper: http://www.larsson.dogma.net/ssrev-tr.pdf
// code:  http://www.larsson.dogma.net/qsufsort.c

// This algorithm computes the suffix array sa by computing its inverse.
// Consecutive groups of suffixes in sa are labeled as sorted groups or
// unsorted groups. For a given pass of the sorter, all suffixes are ordered
// up to their first h characters, and sa is h-ordered. Suffixes in their
// final positions and unambiguously sorted in h-order are in a sorted group.
// Consecutive groups of suffixes with identical first h characters are an
// unsorted group. In each pass of the algorithm, unsorted groups are sorted
// according to the group number of their following suffix.

// In the implementation, if sa[i] is negative, it indicates that i is
// the first element of a sorted group of length -sa[i], and can be skipped.
// An unsorted group sa[i:k] is given the group number of the index of its
// last element, k-1. The group numbers are stored in the inverse slice (inv),
// and when all groups are sorted, this slice is the inverse suffix array.

type Mapper struct {
	sa             []int
	counts         []map[rune]int
	reference      string
	referenceName  string
	firstOccurance map[rune]int
}

func CreateMapper(referenceFile string) Mapper {

	f1, err1 := os.Open(referenceFile)
	if err1 != nil {
		panic(err1)
	}
	defer f1.Close()

	reference := ""

	scanner := bufio.NewScanner(f1)
	scanner.Split(bufio.ScanLines)
	scanner.Scan()
	header := scanner.Text()
	for scanner.Scan() {
		reference += scanner.Text()
	}

	reference += "$"

	reference = strings.ToUpper(reference)

	sa := qsufsort([]byte(reference))

	bwt := make([]rune, len(reference))

	for i, index := range sa {
		if index == 0 {
			bwt[i] = rune(reference[len(reference)-1])
		} else {
			bwt[i] = rune(reference[index-1])
		}
	}

	firstOccurance := make(map[rune]int)
	uniqueCharacters := make([]rune, 0)

	for index, val := range SortString(reference) {
		_, exists := firstOccurance[rune(val)]
		if !exists {
			firstOccurance[rune(val)] = index
			uniqueCharacters = append(uniqueCharacters, rune(val))
		}

	}

	counts := make([]map[rune]int, len(reference)+1)

	for i := 0; i < (len(reference) + 1); i++ {
		counts[i] = make(map[rune]int)
		for _, char := range uniqueCharacters {
			if i == 0 {
				counts[i][char] = 0
			} else {
				counts[i][char] = counts[i-1][char]
				if bwt[i-1] == char {
					counts[i][char]++
				}
			}
		}
	}
	return Mapper{sa, counts, reference, header[1:], firstOccurance}
}

func (m *Mapper) findOccurances(search string) []int {
	top := 0
	bottom := len(m.sa) - 1

	for i := (len(search) - 1); i >= 0; i-- {
		char := rune(search[i])
		top = m.firstOccurance[char] + m.counts[top][char]
		bottom = m.firstOccurance[char] + m.counts[bottom+1][char] - 1
		if top > bottom {
			break
		}
	}
	indices := make([]int, 0)
	if top <= bottom {
		for i := 0; i <= (bottom - top); i++ {
			indices = append(indices, m.sa[top+i])
		}
	}
	return indices
}

func qsufsort(data []byte) []int {
	// initial sorting by first byte of suffix
	sa := sortedByFirstByte(data)
	if len(sa) < 2 {
		return sa
	}
	// initialize the group lookup table
	// this becomes the inverse of the suffix array when all groups are sorted
	inv := initGroups(sa, data)

	// the index starts 1-ordered
	sufSortable := &suffixSortable{sa: sa, inv: inv, h: 1}

	for sa[0] > -len(sa) { // until all suffixes are one big sorted group
		// The suffixes are h-ordered, make them 2*h-ordered
		pi := 0 // pi is first position of first group
		sl := 0 // sl is negated length of sorted groups
		for pi < len(sa) {
			if s := sa[pi]; s < 0 { // if pi starts sorted group
				pi -= s // skip over sorted group
				sl += s // add negated length to sl
			} else { // if pi starts unsorted group
				if sl != 0 {
					sa[pi+sl] = sl // combine sorted groups before pi
					sl = 0
				}
				pk := inv[s] + 1 // pk-1 is last position of unsorted group
				sufSortable.sa = sa[pi:pk]
				sort.Sort(sufSortable)
				sufSortable.updateGroups(pi)
				pi = pk // next group
			}
		}
		if sl != 0 { // if the array ends with a sorted group
			sa[pi+sl] = sl // combine sorted groups at end of sa
		}

		sufSortable.h *= 2 // double sorted depth
	}

	for i := range sa { // reconstruct suffix array from inverse
		sa[inv[i]] = i
	}
	return sa
}

func sortedByFirstByte(data []byte) []int {
	// total byte counts
	var count [256]int
	for _, b := range data {
		count[b]++
	}
	// make count[b] equal index of first occurrence of b in sorted array
	sum := 0
	for b := range count {
		count[b], sum = sum, count[b]+sum
	}
	// iterate through bytes, placing index into the correct spot in sa
	sa := make([]int, len(data))
	for i, b := range data {
		sa[count[b]] = i
		count[b]++
	}
	return sa
}

func initGroups(sa []int, data []byte) []int {
	// label contiguous same-letter groups with the same group number
	inv := make([]int, len(data))
	prevGroup := len(sa) - 1
	groupByte := data[sa[prevGroup]]
	for i := len(sa) - 1; i >= 0; i-- {
		if b := data[sa[i]]; b < groupByte {
			if prevGroup == i+1 {
				sa[i+1] = -1
			}
			groupByte = b
			prevGroup = i
		}
		inv[sa[i]] = prevGroup
		if prevGroup == 0 {
			sa[0] = -1
		}
	}
	// Separate out the final suffix to the start of its group.
	// This is necessary to ensure the suffix "a" is before "aba"
	// when using a potentially unstable sort.
	lastByte := data[len(data)-1]
	s := -1
	for i := range sa {
		if sa[i] >= 0 {
			if data[sa[i]] == lastByte && s == -1 {
				s = i
			}
			if sa[i] == len(sa)-1 {
				sa[i], sa[s] = sa[s], sa[i]
				inv[sa[s]] = s
				sa[s] = -1 // mark it as an isolated sorted group
				break
			}
		}
	}
	return inv
}

type suffixSortable struct {
	sa  []int
	inv []int
	h   int
	buf []int // common scratch space
}

func (x *suffixSortable) Len() int           { return len(x.sa) }
func (x *suffixSortable) Less(i, j int) bool { return x.inv[x.sa[i]+x.h] < x.inv[x.sa[j]+x.h] }
func (x *suffixSortable) Swap(i, j int)      { x.sa[i], x.sa[j] = x.sa[j], x.sa[i] }

func (x *suffixSortable) updateGroups(offset int) {
	bounds := x.buf[0:0]
	group := x.inv[x.sa[0]+x.h]
	for i := 1; i < len(x.sa); i++ {
		if g := x.inv[x.sa[i]+x.h]; g > group {
			bounds = append(bounds, i)
			group = g
		}
	}
	bounds = append(bounds, len(x.sa))
	x.buf = bounds

	// update the group numberings after all new groups are determined
	prev := 0
	for _, b := range bounds {
		for i := prev; i < b; i++ {
			x.inv[x.sa[i]] = offset + b - 1
		}
		if b-prev == 1 {
			x.sa[prev] = -1
		}
		prev = b
	}
}

type sortRunes []rune

func (s sortRunes) Less(i, j int) bool {
	return s[i] < s[j]
}

func (s sortRunes) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

func (s sortRunes) Len() int {
	return len(s)
}

func SortString(s string) string {
	r := []rune(s)
	sort.Sort(sortRunes(r))
	return string(r)
}
