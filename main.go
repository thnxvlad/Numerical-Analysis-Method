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
		func(x float64) float64 {
			return math.Abs(2*math.Sin(2*x) - 1)
		},
	}

	var a, b float64 = -2, 2

	nodesAmount := []int{
		5, 10, 15, 20,
	}

	for i, fn := range F {
		for _, n := range nodesAmount {
			fmt.Printf("#%d EquidistantNodes: %d\n", i, n)
			var eqNip NewtonInterpolationPolynomial
			eqNip.New(fn, gridOfEquidistantNodes(a, b, n))
			fmt.Println(eqNip.String())
			fmt.Printf("Deviation: %f\n", eqNip.CalcDeviation(gridOfEquidistantNodes(a, b, 100)))

			fmt.Printf("#%d ChebyshevNodes: %d\n", i, n)
			var chNip NewtonInterpolationPolynomial
			chNip.New(fn, gridOfChebyshevNodes(a, b, n))
			fmt.Println(chNip.String())
			fmt.Printf("Deviation: %f\n", chNip.CalcDeviation(gridOfEquidistantNodes(a, b, 100)))
		}

	}

}
