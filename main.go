package main

import (
	"fmt"
	"math"
)

func main() {
	const epsilon float64 = 1.0 / 1e6
	const maxIterations int = 100

	A := [][]float64{
		{3, 4, -9, 5},
		{-15, -12, 50, -16},
		{-27, -36, 73, 8},
		{9, 12, -10, -16},
	}
	b := []float64{
		-14, 44, 142, -76,
	}
	fmt.Println(gaussMethod(A, b, epsilon))
	F := []func([]float64) float64{
		func(X []float64) float64 {
			return math.Cos(X[0]-X[1]) - X[0]*X[1] + 2
		},
		func(X []float64) float64 {
			return math.Pow(X[0], 2) + X[0]*X[1] - math.Pow(X[1], 2) + 1.25
		},
	}
	W := [][]func([]float64) float64{
		{
			func(X []float64) float64 {
				return -math.Sin(X[0]-X[1]) - X[1]
			},
			func(X []float64) float64 {
				return math.Sin(X[0]-X[1]) - X[0]
			},
		},
		{
			func(X []float64) float64 {
				return 2*X[0] + X[1]
			},
			func(X []float64) float64 {
				return X[0] - 2*X[1]
			},
		},
	}
	X := []float64{
		-1,
		-2.5,
	}
	fmt.Println(newtonMethod(F, W, X, epsilon, maxIterations))
	fmt.Println(broydenMethod(F, W, X, epsilon, maxIterations))
}
