package main

import (
	"math"
)

// vec * num
func productVectorByNum(vec []float64, num float64) []float64 {
	resultVec := make([]float64, len(vec))
	for i, el := range vec {
		resultVec[i] = num * el
	}
	return resultVec
}

// vec1 + vec2
func vectorSum(vec1 []float64, vec2 []float64) []float64 {
	resultVec := make([]float64, len(vec1))
	for i, _ := range resultVec {
		resultVec[i] = vec1[i] + vec2[i]
	}
	return resultVec
}

// vec1 - vec2
func vectorDifference(vec1 []float64, vec2 []float64) []float64 {
	return vectorSum(vec1, productVectorByNum(vec2, -1))
}

func vectorRounding(vec []float64, epsilon float64) []float64 {
	resultVec := make([]float64, len(vec))
	for i, el := range vec {
		resultVec[i] = math.Round(el/epsilon) / math.Round(1/epsilon)
	}
	return resultVec
}

func areVectorsEqual(vec1 []float64, vec2 []float64, epsilon float64) bool {
	if len(vec1) == len(vec2) {
		for i, _ := range vec1 {
			if math.Abs(vec1[i]-vec2[i]) > epsilon {
				return false
			}
		}
	} else {
		return false
	}
	return true
}

func fnVec(fnVec []func([]float64) float64, X []float64) []float64 {
	resultVec := make([]float64, len(fnVec))
	for i, fn := range fnVec {
		resultVec[i] = fn(X)
	}
	return resultVec
}

func jacobiMatrix(W [][]func([]float64) float64, X []float64) [][]float64 {
	resultMatrix := make([][]float64, len(W))
	for i, fnVector := range W {
		resultMatrix[i] = fnVec(fnVector, X)
	}
	return resultMatrix
}

func matrixProduct(A [][]float64, B [][]float64) [][]float64 {
	resultMatrix := make([][]float64, len(A))
	for i, _ := range A {
		resultMatrix[i] = make([]float64, len(B[0]))
		for j, _ := range A[i] {
			for k, _ := range B[j] {
				resultMatrix[i][k] += A[i][j] * B[j][k]
			}
		}
	}
	return resultMatrix
}

func matrixSum(A [][]float64, B [][]float64) [][]float64 {
	resultMatrix := make([][]float64, len(A))
	for i, _ := range A {
		resultMatrix[i] = make([]float64, len(A[0]))
		for j, _ := range A[i] {
			resultMatrix[i][j] = A[i][j] + B[i][j]
		}
	}
	return resultMatrix
}

func matrixTransposition(A [][]float64) [][]float64 {
	resultMatrix := make([][]float64, len(A[0]))
	for i, _ := range resultMatrix {
		resultMatrix[i] = make([]float64, len(A))
		for j, _ := range resultMatrix[i] {
			resultMatrix[i][j] = A[j][i]
		}
	}
	return resultMatrix
}

func wrapVector(vec []float64) [][]float64 {
	return [][]float64{vec}
}

func getVectorById(A [][]float64, id int) []float64 {
	return A[id]
}

func frobeniusVectorNorm(vec []float64) float64 {
	var result float64 = 0
	for _, el := range vec {
		result += el * el
	}
	return math.Sqrt(result)
}

func productMatrixByNum(A [][]float64, num float64) [][]float64 {
	resultMatrix := make([][]float64, len(A))
	for i, vec := range A {
		resultMatrix[i] = productVectorByNum(vec, num)
	}
	return resultMatrix
}
