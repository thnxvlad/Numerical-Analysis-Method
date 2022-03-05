package main

import (
	"fmt"
	"math"
)

var epsilon float64 = 1.0 / 10000

func calcProductVectorByNum(vec []float64, num float64) []float64 {
	resultVec := make([]float64, len(vec))
	for i, el := range vec {
		resultVec[i] = num * el
	}
	return resultVec
}

func subtractVector(vec1 []float64, vec2 []float64) {
	for i, _ := range vec1 {
		vec1[i] -= vec2[i]
		if math.Abs(vec1[i]) < epsilon {
			vec1[i] = 0
		}
	}
}

func gauss(_A [][]float64, _b []float64) []float64 {
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
		for idx, vec := range A[i+1:] {
			if math.Abs(vec[i]) > pivotElementAbs {
				pivotElementAbs = math.Abs(vec[i])
				pivotElementIdx = idx + i + 1
			}
		}
		A[i], A[pivotElementIdx] = A[pivotElementIdx], A[i]
		b[i], b[pivotElementIdx] = b[pivotElementIdx], b[i]
		for j := i + 1; j < len(A); j++ {
			var scalingMultiplier float64 = A[j][i] / A[i][i]
			subtractVector(A[j], calcProductVectorByNum(A[i], scalingMultiplier))
			b[j] -= scalingMultiplier * b[i]
			if math.Abs(b[j]) < epsilon {
				b[j] = 0
			}
		}
	}

	for i := len(A) - 1; i >= 0; i-- {
		for j := i + 1; j < len(A); j++ {
			b[i] -= b[j] * A[i][j]
		}
		b[i] = math.Round(b[i]/A[i][i]/epsilon) * epsilon
	}
	return b
}

func calcFnVec(fnVec []func(float64, float64) float64, x float64, y float64) []float64 {
	resVec := make([]float64, len(fnVec))
	for i, fn := range fnVec {
		resVec[i] = fn(x, y)
	}
	return resVec
}

func calcJacobiMatrix(W [][]func(float64, float64) float64, x float64, y float64) [][]float64 {
	resMatrix := make([][]float64, len(W))
	for i, fnVec := range W {
		resMatrix[i] = calcFnVec(fnVec, x, y)
	}
	return resMatrix
}

func newtonMethod(F []func(float64, float64) float64, W [][]func(float64, float64) float64, x float64, y float64, iterationsNum int) {
	data := gauss(calcJacobiMatrix(W, x, y), calcProductVectorByNum(calcFnVec(F, x, y), -1))
	newX := x + data[0]
	newY := y + data[1]
	if math.Max(math.Abs(data[0]), math.Abs(data[1])) > epsilon {
		newtonMethod(F, W, newX, newY, iterationsNum+1)
	} else {
		fmt.Println("x: ", newX, ", y: ", newY)
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
	fmt.Println(gauss(A, b))
	F := []func(float64, float64) float64{
		func(x float64, y float64) float64 {
			return math.Sin(2*x-y) - 1.2*x - 0.4
		},
		func(x float64, y float64) float64 {
			return 0.8*math.Pow(x, 2) + 1.5*math.Pow(y, 2) - 1
		},
	}
	W := [][]func(float64, float64) float64{
		{
			func(x float64, y float64) float64 {
				return 2*math.Cos(2*x-y) - 1.2
			},
			func(x float64, y float64) float64 {
				return -math.Cos(2*x - y)
			},
		},
		{
			func(x float64, _ float64) float64 {
				return 1.6 * x
			},
			func(_ float64, y float64) float64 {
				return 3 * y
			},
		},
	}
	newtonMethod(F, W, 0.4, -0.75, 1)
}
