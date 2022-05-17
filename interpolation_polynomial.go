package main

import (
	"fmt"
	"math"
)

// NewtonInterpolationPolynomial - NIP
type NewtonInterpolationPolynomial struct {
	f            func(float64) float64
	gridNodes    []float64
	coefficients []float64
}

func (nip *NewtonInterpolationPolynomial) New(f func(float64) float64, nodes []float64) {
	// SplitDifferenceTable
	table := make([][]float64, len(nodes))

	table[0] = make([]float64, len(nodes))
	for i, node := range nodes {
		table[0][i] = f(node)
	}

	for i := 1; i < len(nodes); i++ {
		table[i] = make([]float64, len(nodes)-i)
		for j := range table[i] {
			table[i][j] = (table[i-1][j+1] - table[i-1][j]) / (nodes[j+i] - nodes[j])
		}
	}

	nip.f = f
	nip.gridNodes = nodes
	nip.coefficients = make([]float64, len(table))
	for i, coeffs := range table {
		nip.coefficients[i] = coeffs[0]
	}
}

func (nip *NewtonInterpolationPolynomial) String() string {
	var polynomial string

	for i, c := range nip.coefficients {
		if i != 0 {
			if c >= 0 {
				polynomial += " + "
			} else {
				polynomial += " - "
			}
		}
		polynomial += fmt.Sprintf("%.8f", math.Abs(c))
		for j := 0; j < i; j++ {
			polynomial += " * (x"
			if nip.gridNodes[j] > 0 {
				polynomial += " - "
			} else {
				polynomial += " + "
			}
			polynomial += fmt.Sprintf("%.8f", math.Abs(nip.gridNodes[j])) + ")"
		}
	}

	return polynomial
}

func (nip *NewtonInterpolationPolynomial) Calc(x float64) float64 {
	var result float64 = 0

	for i, c := range nip.coefficients {
		temp := c
		for j := 0; j < i; j++ {
			temp *= x - nip.gridNodes[j]
		}
		result += temp
	}

	return result
}

func (nip *NewtonInterpolationPolynomial) CalcDeviation(nodes []float64) float64 {
	maxDeviation := 0.0

	for _, x := range nodes {
		maxDeviation = math.Max(maxDeviation, math.Abs(nip.Calc(x)-nip.f(x)))
	}

	return maxDeviation
}

func gridOfEquidistantNodes(a, b float64, n int) []float64 {
	grid := make([]float64, n+1)

	for i := range grid {
		grid[i] = a + float64(i)*(b-a)/float64(n)
	}

	return grid
}

func gridOfChebyshevNodes(a, b float64, n int) []float64 {
	grid := make([]float64, n+1)

	for i := range grid {
		grid[i] = (a+b)/2 + (b-a)/2*math.Cos(float64(2*i+1)*math.Pi/float64(2*(n+1)))
	}

	return grid
}
