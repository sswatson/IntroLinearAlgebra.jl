using IntroLinearAlgebra
using Test

# write your own tests here
@test rowswitch([1 2 3 4; 5 6 7 8],1,2) == [5 6 7 8; 1 2 3 4]

@test rowscale([1 2 3; 4 5 6],2,-3//2) == [1 2 3; -6 -15//2 -9]

@test rowadd([0.3 0.6; π 2π],2,1,-π/0.3) == [0.3 0.6; 0.0 0.0]

@test rref([1 2 3 4 5; 6 7 8 9 10]) == [1 0 -1 -2 -3; 0 1 2 3 4]

@test length(rref(randn(8,8))) == 64

A = [1//1 2 3 4;
        5 6 7 8]

zerooutbelow!(A,1,1)

@test A == [1 2 3 4; 0 -4 -8 -12]

