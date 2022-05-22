package main

import (
	"fmt"
	"math"
)

type LeastSquaresApproximationPolynomial struct {
	f            func(float64) float64
	coefficients []float64
}

func (p *LeastSquaresApproximationPolynomial) New(f func(float64) float64, nodes []float64, n int, epsilon float64) error {
	s := make([]float64, 2*n+1)
	for i := range s {
		s[i] = 0
		for _, x := range nodes {
			s[i] += math.Pow(x, float64(i))
		}
	}

	m := make([]float64, n+1)
	for i := range m {
		m[i] = 0
		for _, x := range nodes {
			m[i] += f(x) * math.Pow(x, float64(i))
		}
	}

	matrix := make([][]float64, n+1)
	for i := 0; i <= n; i++ {
		matrix[i] = make([]float64, n+1)
		for j := 0; j <= n; j++ {
			matrix[i][j] = s[j+i]
		}
	}

	c, err := gaussMethod(matrix, m, epsilon)
	if err != nil {
		return err
	}

	p.coefficients = c
	p.f = f

	return nil
}

func (p *LeastSquaresApproximationPolynomial) Calc(x float64) float64 {
	result := 0.0
	for i, c := range p.coefficients {
		result += c * math.Pow(x, float64(i))
	}

	return result
}

func (p *LeastSquaresApproximationPolynomial) String() string {
	var polynomial string
	for i, c := range p.coefficients {
		if i == 0 {
			polynomial = fmt.Sprintf("%.3f", c)
		}

		if c >= 0 {
			polynomial += " + "
		} else {
			polynomial += " - "
		}

		polynomial += fmt.Sprintf("%.3f * x^%d", math.Abs(c), i)
	}

	return polynomial
}

func (p *LeastSquaresApproximationPolynomial) CalcDeviation(nodes []float64) float64 {
	dev := 0.0

	for _, x := range nodes {
		dev += math.Pow(p.f(x)-p.Calc(x), 2)
	}

	return dev
}
