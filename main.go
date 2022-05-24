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

	dF := []func(float64) float64{
		func(x float64) float64 {
			return math.Exp(math.Sin(x)) * math.Cos(x)
		},
	}

	var a, b, epsilon float64 = -2, 2, 0.0000000001

	N := []int{
		10,
	}

	for i, fn := range F {
		for j := range N {
			fmt.Printf("#%d n: %d\n", i, N[j])
			var s CubicSpline
			err := s.New(fn, dF[i], gridOfEquidistantNodes(a, b, N[j]), epsilon)
			if err != nil {
				fmt.Println(err)
				continue
			}
			dev, err := s.CalcDeviation(gridOfEquidistantNodes(a, b, 100))
			if err != nil {
				fmt.Println(err)
				continue
			}
			fmt.Printf("Deviation: %f\n", dev)
		}
	}
}
