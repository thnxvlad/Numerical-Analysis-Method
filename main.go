package main

import (
	"errors"
	"fmt"
	"math"
)

const EPSILON float64 = 1.0 / 10000
const MAX_ITERATIONS int = 100

// vec * num
func calcProductVectorByNum(vec []float64, num float64) []float64 {
	resultVec := make([]float64, len(vec))
	for i, el := range vec {
		resultVec[i] = num * el
	}
	return resultVec
}

// vec1 + vec2
func calcVectorSum(vec1 []float64, vec2 []float64) []float64 {
	resultVec := make([]float64, len(vec1))
	for i, _ := range resultVec {
		resultVec[i] = vec1[i] + vec2[i]
	}
	return resultVec
}

// vec1 - vec2
func calcVectorDifference(vec1 []float64, vec2 []float64) []float64 {
	return calcVectorSum(vec1, calcProductVectorByNum(vec2, -1))
}

func calcVectorRounding(vec []float64) []float64 {
	resultVec := make([]float64, len(vec))
	for i, el := range vec {
		resultVec[i] = math.Round(el/EPSILON) * EPSILON
	}
	return resultVec
}

// Ax = b
func gaussMethod(_A [][]float64, _b []float64) []float64 {
	A := make([][]float64, len(_A))
	for i, slice := range _A {
		A[i] = make([]float64, len(slice))
		copy(A[i], slice)
	}
	b := make([]float64, len(_b))
	copy(b, _b)
	for i := 0; i < len(A); i++ {
		var pivotElementAbs float64 = math.Abs(A[i][i])
		var pivotElementIdx int = i
		for j, vec := range A[i+1:] {
			if math.Abs(vec[i]) > pivotElementAbs {
				pivotElementAbs = math.Abs(vec[i])
				pivotElementIdx = j + i + 1
			}
		}
		A[i], A[pivotElementIdx] = A[pivotElementIdx], A[i]
		b[i], b[pivotElementIdx] = b[pivotElementIdx], b[i]
		for j := i + 1; j < len(A); j++ {
			var scalingMultiplier float64 = A[j][i] / A[i][i]
			A[j] = calcVectorDifference(A[j], calcProductVectorByNum(A[i], scalingMultiplier))
			b[j] -= scalingMultiplier * b[i]
		}
	}

	for i := len(A) - 1; i >= 0; i-- {
		for j := i + 1; j < len(A); j++ {
			b[i] -= b[j] * A[i][j]
		}
		b[i] /= A[i][i]
	}
	b = calcVectorRounding(b)
	return b
}

func calcFnVec(fnVec []func([]float64) float64, X []float64) []float64 {
	resultVec := make([]float64, len(fnVec))
	for i, fn := range fnVec {
		resultVec[i] = fn(X)
	}
	return resultVec
}

func calcJacobiMatrix(W [][]func([]float64) float64, X []float64) [][]float64 {
	resultMatrix := make([][]float64, len(W))
	for i, fnVec := range W {
		resultMatrix[i] = calcFnVec(fnVec, X)
	}
	return resultMatrix
}

func newtonMethod(F []func([]float64) float64, W [][]func([]float64) float64, X []float64) ([]float64, int, error) {
	newX := X
	var iteration int = 0
	for {
		if iteration == MAX_ITERATIONS {
			return calcVectorRounding(newX), iteration, errors.New("iteration limit exceeded")
		} else {
			iteration++
		}
		dX := gaussMethod(calcJacobiMatrix(W, newX), calcProductVectorByNum(calcFnVec(F, newX), -1))
		newX = calcVectorSum(newX, dX)
		var calcContinueCondition bool = false
		for _, dx := range dX {
			if dx > EPSILON {
				calcContinueCondition = true
				break
			}
		}
		switch calcContinueCondition {
		case true:
			continue
		case false:
			return calcVectorRounding(newX), iteration, nil
		}
	}
}

func main() {
	A := [][]float64{
		{3, 4, -9, 5},
		{-15, -12, 50, -16},
		{-27, -36, 73, 8},
		{9, 12, -10, -16},
	}
	b := []float64{
		-14, 44, 142, -76,
	}
	fmt.Println(gaussMethod(A, b))
	F := []func([]float64) float64{
		func(X []float64) float64 {
			return math.Sin(2*X[0]-X[1]) - 1.2*X[0] - 0.4
		},
		func(X []float64) float64 {
			return 0.8*math.Pow(X[0], 2) + 1.5*math.Pow(X[1], 2) - 1
		},
	}
	W := [][]func([]float64) float64{
		{
			func(X []float64) float64 {
				return 2*math.Cos(2*X[0]-X[1]) - 1.2
			},
			func(X []float64) float64 {
				return -math.Cos(2*X[0] - X[1])
			},
		},
		{
			func(X []float64) float64 {
				return 1.6 * X[0]
			},
			func(X []float64) float64 {
				return 3 * X[1]
			},
		},
	}
	X := []float64{
		0.4,
		-0.75,
	}
	fmt.Println(newtonMethod(F, W, X))
}
