package main

import (
	"errors"
	"math"
)

// Ax = b
func gaussMethod(_A [][]float64, _b []float64, epsilon float64) ([]float64, error) {
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
			A[j] = vectorDifference(A[j], productVectorByNum(A[i], scalingMultiplier))
			b[j] -= scalingMultiplier * b[i]
		}
	}

	for i := len(A) - 1; i >= 0; i-- {
		for j := i + 1; j < len(A); j++ {
			b[i] -= b[j] * A[i][j]
		}
		if math.Abs(A[i][i]) < epsilon {
			return vectorRounding(b, epsilon), errors.New("the vector system is linearly dependent")
		}
		b[i] /= A[i][i]
	}
	b = vectorRounding(b, epsilon)
	return b, nil
}

func newtonMethod(F []func([]float64) float64, W [][]func([]float64) float64, X []float64, epsilon float64, maxIterations int) ([]float64, error) {
	newX := X
	var iteration int = 0
	for {
		if iteration == maxIterations {
			return vectorRounding(newX, epsilon), errors.New("iteration limit exceeded")
		} else {
			iteration++
		}
		dX, err := gaussMethod(jacobiMatrix(W, newX), productVectorByNum(fnVec(F, newX), -1), epsilon)
		if err != nil {
			return vectorRounding(newX, epsilon), err
		}
		newX = vectorSum(newX, dX)
		var continueCondition bool = false
		for _, dx := range dX {
			if math.Abs(dx) > epsilon {
				continueCondition = true
				break
			}
		}
		switch continueCondition {
		case true:
			continue
		case false:
			return vectorRounding(newX, epsilon), nil
		}
	}
}

func broydenMethod(F []func([]float64) float64, W [][]func([]float64) float64, X []float64, epsilon float64, maxIterations int) ([]float64, error) {
	J := jacobiMatrix(W, X)
	newX := X
	oldX := X
	var iteration int = 0
	for {
		if iteration == maxIterations {
			return vectorRounding(newX, epsilon), errors.New("iteration limit exceeded")
		} else {
			iteration++
		}
		dX, err := gaussMethod(J, productVectorByNum(fnVec(F, newX), -1), epsilon)
		if err != nil {
			return vectorRounding(newX, epsilon), err
		}
		oldX = newX
		newX = vectorSum(newX, dX)
		dF := vectorDifference(fnVec(F, newX), fnVec(F, oldX))
		J = matrixSum(
			J,
			productMatrixByNum(
				matrixProduct(
					matrixTransposition(wrapVector(vectorDifference(
						dF,
						getVectorById(matrixTransposition(matrixProduct(
							J,
							matrixTransposition(wrapVector(dX)),
						)), 0),
					))),
					wrapVector(dX),
				),
				1.0/math.Pow(frobeniusVectorNorm(dX), 2),
			),
		)
		var continueCondition bool = false
		for _, dx := range dX {
			if math.Abs(dx) > epsilon {
				continueCondition = true
				break
			}
		}
		switch continueCondition {
		case true:
			continue
		case false:
			return vectorRounding(newX, epsilon), nil
		}
	}
}
