package main

import (
	"fmt"
	"math"
)

func main() {
	F := []func(float64) float64{
		func(x float64) float64 {
			return math.Exp(math.Sin(x))
		},
	}

	var a, b, h, epsilon float64 = -2, 2, 0.1, 0.0000001

	n := []int{
		3, 5,
	}

	for i, fn := range F {
		for j := range n {
			fmt.Printf("#%d n: %d\n", i, n[j])
			N := int((b - a) / h)
			var p LeastSquaresApproximationPolynomial
			err := p.New(fn, gridOfEquidistantNodes(a, b, N), n[j], epsilon)
			if err != nil {
				fmt.Println(err)
				continue
			}
			fmt.Println(p.String())
			fmt.Printf("Deviation: %f\n", p.CalcDeviation(gridOfEquidistantNodes(a, b, N)))
		}
	}
}
