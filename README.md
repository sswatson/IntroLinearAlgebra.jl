# IntroLinearAlgebra

`IntroLinearAlgebra` is a Julia package designed for students in an introductory linear algebra course. It uses MathJax to display matrices in an IJulia notebook, and it provides functions for step-by-step application of elementary row operations.

## Installation

Run `Pkg.clone("git://github.com/sswatson/IntroLinearAlgebra.jl.git")` from a Julia prompt. 

## Examples

```julia
julia> using IntroLinearAlgebra
julia> M = [1 2 3; 4 5 6; 7 8 9]
3×3 Array{Int64,2}:
 1  2  3
 4  5  6
 7  8  9

julia> rref(M)
3×3 Array{Rational{Int64},2}:
 1//1  0//1  -1//1
 0//1  1//1   2//1
 0//1  0//1   0//1

julia> M = rowadd(M,2,1,-4)
3×3 Array{Int64,2}:
 1   2   3
 0  -3  -6
 7   8   9

julia> M = rowswitch(M,3,1)
3×3 Array{Int64,2}:
 7   8   9
 0  -3  -6
 1   2   3

julia> M = rowadd(M,3,1,-1//7)
3×3 Array{Rational{Int64},2}:
 7//1   8//1   9//1
 0//1  -3//1  -6//1
 0//1   6//7  12//7

julia> rowscale(M,1,1//7)
3×3 Array{Rational{Int64},2}:
 1//1   8//7   9//7
 0//1  -3//1  -6//1
 0//1   6//7  12//7

julia> setdigits(8) #changes the default number of digits to display in IJulia
```

[![Build Status](https://travis-ci.org/sswatson/IntroLinearAlgebra.jl.svg?branch=master)](https://travis-ci.org/sswatson/IntroLinearAlgebra.jl)

[![Coverage Status](https://coveralls.io/repos/sswatson/IntroLinearAlgebra.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/sswatson/IntroLinearAlgebra.jl?branch=master)

[![codecov.io](http://codecov.io/github/sswatson/IntroLinearAlgebra.jl/coverage.svg?branch=master)](http://codecov.io/github/sswatson/IntroLinearAlgebra.jl?branch=master)
