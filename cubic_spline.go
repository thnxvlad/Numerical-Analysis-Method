package main

import (
	"errors"
	"math"
)

type CubicSpline struct {
	f            func(float64) float64
	coefficients []float64
	nodes        []float64
}

func (s *CubicSpline) buildMVector(f func(float64) float64, df func(float64) float64, nodes []float64, epsilon float64) ([]float64, error) {
	N := len(nodes) - 1

	A := make([][]float64, N+1)
	for i := range A {
		A[i] = make([]float64, N+1)
	}

	fM := make([]float64, N+1)

	h := nodes[1] - nodes[0]

	A[0][0] = h / 3
	A[0][1] = h / 6

	fM[0] = (f(nodes[1])-f(nodes[0]))/h - df(nodes[0])

	for i := 1; i < N; i++ {
		A[i][i-1] = h / 6
		A[i][i] = 2 * h / 3
		A[i][i+1] = h / 6
		fM[i] = (f(nodes[i+1])-f(nodes[i]))/h - (f(nodes[i])-f(nodes[i-1]))/h
	}

	A[N][N-1] = h / 6
	A[N][N] = h / 3
	fM[N] = df(nodes[N]) - (f(nodes[N])-f(nodes[N-1]))/h

	M, err := tridiagonalMatrixAlgorithm(A, fM, epsilon)
	if err != nil {
		return nil, err
	}

	return M, nil
}

func (s *CubicSpline) New(f func(float64) float64, df func(float64) float64, nodes []float64, epsilon float64) error {
	M, err := s.buildMVector(f, df, nodes, epsilon)
	if err != nil {
		return err
	}

	s.f = f
	s.coefficients = M
	s.nodes = nodes

	return nil
}

func (s *CubicSpline) Calc(x float64) (float64, error) {
	h := s.nodes[1] - s.nodes[0]
	if x < s.nodes[0] {
		return math.MaxFloat64, errors.New("x value less than bounds of spline")
	}
	if x > s.nodes[len(s.nodes)-1] {
		return math.MaxFloat64, errors.New("x value more than bounds of spline")
	}

	j := int((x - s.nodes[0]) / h)
	if x != s.nodes[len(s.nodes)-1] {
		j++
	}
	res := s.coefficients[j-1]*math.Pow(s.nodes[j]-x, 3)/(6*h) +
		s.coefficients[j]*math.Pow(x-s.nodes[j-1], 3)/(6*h) +
		(s.f(s.nodes[j-1])-s.coefficients[j-1]*math.Pow(h, 2)/6)*(s.nodes[j]-x)/h +
		(s.f(s.nodes[j])-s.coefficients[j]*math.Pow(h, 2)/6)*(x-s.nodes[j-1])/h

	return res, nil
}

func (s *CubicSpline) CalcDeviation(nodes []float64) (float64, error) {
	var maxDeviation float64 = -1

	for _, x := range nodes {
		S, err := s.Calc(x)
		if err != nil {
			return -1, err
		}
		maxDeviation = math.Max(maxDeviation, math.Abs(S-s.f(x)))
	}

	return maxDeviation, nil
}
