package main

import (
	"fmt"
	"math"
	"strconv"
	"testing"
)

func TestGaussMethodTabelDriven(t *testing.T) {
	var tests = []struct {
		A       [][]float64
		b       []float64
		epsilon float64
		x       []float64
	}{
		{
			[][]float64{
				{3, 4, -9, 5},
				{-15, -12, 50, -16},
				{-27, -36, 73, 8},
				{9, 12, -10, -16},
			},
			[]float64{
				-14, 44, 142, -76,
			},
			1.0 / 10000,
			[]float64{
				-8, -2, -2, 0,
			},
		},
		{
			[][]float64{
				{1, -2, 3},
				{3, 1, -1},
				{2, 5, 2},
			},
			[]float64{
				2, 3, 9,
			},
			1.0 / 10000,
			[]float64{
				1, 1, 1,
			},
		},
		{
			[][]float64{
				{1, -1},
				{2, 1},
			},
			[]float64{
				-5, -7,
			},
			1.0 / 10000,
			[]float64{
				-4, 1,
			},
		},
		{
			[][]float64{
				{2, 5, 4, 1},
				{1, 3, 2, 1},
				{2, 10, 9, 7},
				{3, 8, 9, 2},
			},
			[]float64{
				20, 11, 40, 37,
			},
			1.0 / 100,
			[]float64{
				1, 2, 2, 0,
			},
		},
		{
			[][]float64{
				{2},
			},
			[]float64{
				-4,
			},
			1.0 / 10,
			[]float64{
				-2,
			},
		},
	}
	for testIdx, tt := range tests {
		t.Run(strconv.Itoa(testIdx), func(t *testing.T) {
			ans, err := gaussMethod(tt.A, tt.b, tt.epsilon)
			if !areVectorsEqual(ans, tt.x, tt.epsilon) {
				t.Errorf("got: %s expected: %s error: %s", fmt.Sprintln(ans), fmt.Sprintln(tt.x), fmt.Sprintln(err))
			}
		})
	}
}

func TestNewtonMethodTabelDriven(t *testing.T) {
	var tests = []struct {
		F             []func([]float64) float64
		W             [][]func([]float64) float64
		X             []float64
		epsilon       float64
		maxIterations int
		expectedX     []float64
	}{
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Sin(2*X[0]-X[1]) - 1.2*X[0] - 0.4
				},
				func(X []float64) float64 {
					return 0.8*math.Pow(X[0], 2) + 1.5*math.Pow(X[1], 2) - 1
				},
			},
			[][]func([]float64) float64{
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
			},
			[]float64{
				0.4,
				-0.75,
			},
			1.0 / 1e4,
			100,
			[]float64{
				0.4912,
				-0.7335,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) - 1
				},
				func(X []float64) float64 {
					return math.Pow(X[0], 3) - X[1]
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
				},
				{
					func(X []float64) float64 {
						return 3 * math.Pow(X[0], 2)
					},
					func(X []float64) float64 {
						return -1
					},
				},
			},
			[]float64{
				0.9,
				0.5,
			},
			1.0 / 1e4,
			100,
			[]float64{
				0.826,
				0.5636,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return X[0] + X[1] - 3
				},
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) - 9
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 1
					},
					func(X []float64) float64 {
						return 1
					},
				},
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
				},
			},
			[]float64{
				1,
				5,
			},
			1.0 / 1e3,
			100,
			[]float64{
				0,
				3,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) + math.Pow(X[2], 2) - 1
				},
				func(X []float64) float64 {
					return 2*math.Pow(X[0], 2) + math.Pow(X[1], 2) - 4*X[2]
				},
				func(X []float64) float64 {
					return 3*math.Pow(X[0], 2) - 4*X[1] + math.Pow(X[2], 2)
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
					func(X []float64) float64 {
						return 2 * X[2]
					},
				},
				{
					func(X []float64) float64 {
						return 4 * X[0]
					},
					func(X []float64) float64 {
						return 4 * X[1]
					},
					func(X []float64) float64 {
						return -4
					},
				},
				{
					func(X []float64) float64 {
						return 6 * X[0]
					},
					func(X []float64) float64 {
						return -4
					},
					func(X []float64) float64 {
						return 2 * X[2]
					},
				},
			},
			[]float64{
				0.5,
				0.5,
				0.5,
			},
			5.0 / 1e3,
			100,
			[]float64{
				0.785,
				0.495,
				0.37,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return X[0] + 3*math.Log10(X[0]) - math.Pow(X[1], 2)
				},
				func(X []float64) float64 {
					return 2*math.Pow(X[0], 2) - X[0]*X[1] - 5*X[0] + 1
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 1 + 3*0.43429/X[0] - 2*X[1]
					},
					func(X []float64) float64 {
						return -2 * X[1]
					},
				},
				{
					func(X []float64) float64 {
						return 4*X[0] - X[1] - 5
					},
					func(X []float64) float64 {
						return -X[1]
					},
				},
			},
			[]float64{
				3.5,
				2.2,
			},
			1.0 / 1e5,
			100,
			[]float64{
				3.48744,
				2.26162,
			},
		},
	}
	for testIdx, tt := range tests {
		t.Run(strconv.Itoa(testIdx), func(t *testing.T) {
			ans, err := newtonMethod(tt.F, tt.W, tt.X, tt.epsilon, tt.maxIterations)
			if !areVectorsEqual(ans, tt.expectedX, tt.epsilon) {
				t.Errorf("got: %s expected: %s error: %s", fmt.Sprintln(ans), fmt.Sprintln(tt.expectedX), fmt.Sprintln(err))
			}
		})
	}
}

func TestBroydenMethodTabelDriven(t *testing.T) {
	var tests = []struct {
		F             []func([]float64) float64
		W             [][]func([]float64) float64
		X             []float64
		epsilon       float64
		maxIterations int
		expectedX     []float64
	}{
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Sin(2*X[0]-X[1]) - 1.2*X[0] - 0.4
				},
				func(X []float64) float64 {
					return 0.8*math.Pow(X[0], 2) + 1.5*math.Pow(X[1], 2) - 1
				},
			},
			[][]func([]float64) float64{
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
			},
			[]float64{
				0.4,
				-0.75,
			},
			1.0 / 1e4,
			100,
			[]float64{
				0.4912,
				-0.7335,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) - 1
				},
				func(X []float64) float64 {
					return math.Pow(X[0], 3) - X[1]
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
				},
				{
					func(X []float64) float64 {
						return 3 * math.Pow(X[0], 2)
					},
					func(X []float64) float64 {
						return -1
					},
				},
			},
			[]float64{
				0.9,
				0.5,
			},
			1.0 / 1e4,
			100,
			[]float64{
				0.826,
				0.5636,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return X[0] + X[1] - 3
				},
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) - 9
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 1
					},
					func(X []float64) float64 {
						return 1
					},
				},
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
				},
			},
			[]float64{
				1,
				5,
			},
			1.0 / 1e3,
			100,
			[]float64{
				0,
				3,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return math.Pow(X[0], 2) + math.Pow(X[1], 2) + math.Pow(X[2], 2) - 1
				},
				func(X []float64) float64 {
					return 2*math.Pow(X[0], 2) + math.Pow(X[1], 2) - 4*X[2]
				},
				func(X []float64) float64 {
					return 3*math.Pow(X[0], 2) - 4*X[1] + math.Pow(X[2], 2)
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 2 * X[0]
					},
					func(X []float64) float64 {
						return 2 * X[1]
					},
					func(X []float64) float64 {
						return 2 * X[2]
					},
				},
				{
					func(X []float64) float64 {
						return 4 * X[0]
					},
					func(X []float64) float64 {
						return 4 * X[1]
					},
					func(X []float64) float64 {
						return -4
					},
				},
				{
					func(X []float64) float64 {
						return 6 * X[0]
					},
					func(X []float64) float64 {
						return -4
					},
					func(X []float64) float64 {
						return 2 * X[2]
					},
				},
			},
			[]float64{
				0.5,
				0.5,
				0.5,
			},
			5.0 / 1e3,
			100,
			[]float64{
				0.785,
				0.495,
				0.37,
			},
		},
		{
			[]func([]float64) float64{
				func(X []float64) float64 {
					return X[0] + 3*math.Log10(X[0]) - math.Pow(X[1], 2)
				},
				func(X []float64) float64 {
					return 2*math.Pow(X[0], 2) - X[0]*X[1] - 5*X[0] + 1
				},
			},
			[][]func([]float64) float64{
				{
					func(X []float64) float64 {
						return 1 + 3*0.43429/X[0] - 2*X[1]
					},
					func(X []float64) float64 {
						return -2 * X[1]
					},
				},
				{
					func(X []float64) float64 {
						return 4*X[0] - X[1] - 5
					},
					func(X []float64) float64 {
						return -X[1]
					},
				},
			},
			[]float64{
				3.5,
				2.2,
			},
			1.0 / 1e5,
			100,
			[]float64{
				3.48744,
				2.26162,
			},
		},
	}
	for testIdx, tt := range tests {
		t.Run(strconv.Itoa(testIdx), func(t *testing.T) {
			ans, err := broydenMethod(tt.F, tt.W, tt.X, tt.epsilon, tt.maxIterations)
			if !areVectorsEqual(ans, tt.expectedX, tt.epsilon) {
				t.Errorf("got: %s expected: %s error: %s", fmt.Sprintln(ans), fmt.Sprintln(tt.expectedX), fmt.Sprintln(err))
			}
		})
	}
}
