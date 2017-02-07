module IntroLinearAlgebra

# package code goes here

export rowswitch,
       rowscale,
       rowadd,
       rowswitch!,
       rowscale!,
       rowadd!,
       rref,
       zerooutabove!,
       zerooutbelow!,
       normalizerow!,
       setdigits,
       setalignment

import Base.show

if Pkg.installed("SymPy") != nothing
    import SymPy
    sympyexists = true
else 
    sympyexists = false
end

RatOrInt = Union{Rational,Integer}
RatOrIntOrSym = sympyexists ? Union{Rational,Integer,SymPy.Sym} : RatOrInt
RealOrSym = sympyexists ? Union{Real,SymPy.Sym} : Real

DIGITS = 4
ALIGNMENT = true

"""
    setdigits(n)

Set the number of digits after the decimal point
to display, for floating-point matrices
"""
function setdigits(n::Integer)
    global DIGITS
    DIGITS = n
end

"""
    setalignment(b::Bool)

Determine whether to align displayed matrix entries
disregarding the sign (true) or including the sign (false)
"""
function setalignment(b::Bool)
    global ALIGNMENT
    ALIGNMENT = b
end

"""
    texstring(M)

Return a string for processing by LaTeX to pretty print
the matrix `M`
"""
function texstring{T<:Real}(M::Array{T,2})
    ϵ = T<:RatOrInt ? 0 : eps(norm(M,Inf))
    function neg(x::Real) # test whether a number is
                          # genuinely negative (-0.0 isn't)
        return x < -ϵ
    end
    function prettystring(x::RatOrInt;padding=true)
        if isa(x,Integer)
            if padding
                return neg(x) ? string(x) : "\\hphantom{-}"*string(x)
            else
                return string(x)
            end
        elseif x.den == 1
            return prettystring(x.num,padding=padding)
        else
            sgn = neg(x.num) ? "-" : (padding?"\\hphantom{-}":"")
            return sgn*"\\frac{"*string(abs(x.num))*"}{"*string(x.den)*"}"
        end
    end
    function prettystring(x::Real;padding=true)
        global DIGITS
        return (neg(x) ? "-" : (padding ? "\\hphantom{-}" : "")) *
            string(round(abs(x),DIGITS))
    end
    s = "\$\\left[\\begin{array}{" * repeat("c",size(M,2)) * "}"
    global ALIGNMENT
    columnpadding = any(M .< -ϵ, 1) # checks for negatives in each column
    for i=1:size(M,1)
        for j=1:size(M,2)
            thiscolumnpadding = any(M[:,j] .< -ϵ)
            s *= prettystring(M[i,j],padding=ALIGNMENT && columnpadding[j])
            if j < size(M,2)
                s *= " & "
            end
        end
        s *= " \\\\ "
    end
    s *= "\\end{array}\\right]\$"
    return s
end

show{T<:Real}(io::IO, ::MIME"text/latex", s::Array{T,2}) = write(io, texstring(s))

"""
    rowswitch(M,i,j)

Return matrix obtained by switching rows `i` and `j` in matrix `M`

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> rowswitch(M,1,2)
2×3 Array{Rational{Int64},2}:
 4//1  5//1  6//1
 1//1  2//1  3//1
```
"""
function rowswitch{T<:Real}(M::Array{T,2},i::Integer,j::Integer)
    A = copy(M)
    A[i,:], A[j,:] = A[j,:], A[i,:] 
    return A
end

"""
    rowscale(M,i,x)

Return matrix obtained by scaling row `i` in matrix `M` by 
a factor of `x`

```julia 
julia> M = [1 2 3; 4 5 6]
julia> rowscale(M,1,16)
2×3 Array{Int64,2}:
 16  32  48
  4   5   6
```
"""
function rowscale{T<:RatOrIntOrSym}(M::Array{T,2},i::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
    A[i,:] *= x
    return A
end

function rowscale{T<:Real}(M::Array{T,2},i::Integer,x::Real)
    A = copy(M)
    A[i,:] *= x
    return A
end

"""
    rowadd(M,i,j,x)

Return matrix obtained by adding `x` times row `j` to 
row `i` in matrix `M`

```julia 
julia> M = [1 2 3; 4 5 6]
julia> rowadd(M,1,2,-1)
2×3 Array{Int64,2}:
 -3  -3  -3
  4   5   6
```
"""
function rowadd{T<:RatOrIntOrSym}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
    A[i,:] += x*A[j,:]
    return A
end

function rowadd{T<:Real}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    A = copy(M)
    A[i,:] += x*A[j,:]
    return A
end

"""
    rowswitch!(M,i,j)

Version of `rowswitch` that modifies the argument rather than 
returning a new matrix
"""
function rowswitch!{T<:Real}(M::Array{T,2},i::Integer,j::Integer)
    if i == j
        return nothing
    else
        M[i,:], M[j,:] = M[j,:], M[i,:]
        return nothing
    end
end

"""
    rowscale!(M,i,j)

Version of `rowscale` that modifies the argument rather than 
returning a new matrix
"""
function rowscale!{T<:RatOrIntOrSym}(M::Array{T,2},i::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
                    "to use rowscale! with a rational, like M = [1//1 2; 3 4]")
    end
    M[i,:] *= x
    return nothing
end

function rowscale!{T<:Real}(M::Array{T,2},i::Integer,x::Real)
    M[i,:] *= x
    return nothing
end

"""
    rowadd!(M,i,j)

Version of `rowadd` that modifies the argument rather than 
returning a new matrix
"""
function rowadd!{T<:RatOrIntOrSym}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
                    "to use rowadd! with a rational, like M = [1//1 2; 3 4]")
    else
        M[i,:] += x*M[j,:]
    end
    return nothing
end

function rowadd!{T<:Real}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    M[i,:] += x*M[j,:]
    return nothing
end

"""
    rref(M;showsteps=false)

Compute the reduced row echelon form of M and return 
the RREF if `showsteps` is `false`, and the sequence of 
steps performed to compute the RREF if `showsteps` is 
`true`

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> rref(M)
2×3 Array{Rational{Int64},2}:
 1//1  0//1  -1//1
 0//1  1//1   2//1
julia> rref(M;showsteps=true)
4-element Array{Array{Rational{Int64},2},1}:
 Rational{Int64}[1//1 2//1 3//1; 4//1 5//1 6//1]  
 Rational{Int64}[1//1 2//1 3//1; 0//1 -3//1 -6//1]
 Rational{Int64}[1//1 2//1 3//1; 0//1 1//1 2//1]  
 Rational{Int64}[1//1 0//1 -1//1; 0//1 1//1 2//1] 
```
"""
function rref{T<:Real}(M::Array{T,2};showsteps=false)
    A = T<: RatOrInt ? copy(convert(Array{Rational{Int64},2},M)) : copy(M)
    ϵ = T <: RatOrInt ? 0 : eps(norm(A,Inf)) 
    steps = typeof(A)[]
    current_row = 1
    for j=1:size(A,2)
        if all(abs(A[current_row:end,j]) .≤ ϵ) 
            if ϵ > 0
                A[current_row:end,j] = 0
            end
        else
            i = current_row-1 + findfirst(abs(A[current_row:end,j]) .> ϵ)
            rowswitch!(A,current_row,i)
            if showsteps push!(steps,copy(A)) end
            zerooutbelow!(A,current_row,j)
            if showsteps push!(steps,copy(A)) end
            normalizerow!(A,current_row)
            if showsteps push!(steps,copy(A)) end
            current_row += 1
            if current_row > size(A,1)
                break
            end
        end
    end
    for i=size(A,1):-1:1
        if all(abs(A[i,:]) .≤ ϵ)
            continue
        else
            j = findfirst(abs(A[i,:]) .> ϵ) 
            zerooutabove!(A,i,j)
            if showsteps push!(steps,copy(A)) end
        end
    end
    if showsteps
        S = T <: RatOrInt ? Rational{Int64} : T
        uniquesteps = Array{T,2}[]
        for i=1:length(steps)
            if i == 1 || steps[i] != uniquesteps[end]
                push!(uniquesteps,steps[i])
            end
        end
        return uniquesteps
    else
        return A
    end
end

"""
    zerooutbelow!(M,i,j)

Use row operations to zero out entries in column `j` below row `i`

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> zerooutbelow!(M,1,2)
julia> M
2×3 Array{Rational{Int64},2}:
 1//1  2//1   3//1
 3//2  0//1  -3//2
```
"""
function zerooutbelow!{T<:RealOrSym}(M::Array{T,2},i::Integer,j::Integer) # zero out entries below the (i,j)th
    for k=i+1:size(M,1)
        rowadd!(M,k,i,T<:RatOrIntOrSym ? -M[k,j]//M[i,j] : -M[k,j]/M[i,j])
    end
    return nothing
end

#----------------------------------------------------------------
"""
    zerooutabove!(M,i,j)

Use row operations to zero out entries in column `j` above row `i`

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> zerooutabove!(M,2,2)
julia> M
julia> M
2×3 Array{Rational{Int64},2}:
 -3//5  0//1  3//5
  4//1  5//1  6//1
```
"""
function zerooutabove!{T<:RealOrSym}(M::Array{T,2},i::Integer,j::Integer)::Void
    for k=1:i-1
        rowadd!(M,k,i,T<:RatOrIntOrSym ? -M[k,j]//M[i,j] : -M[k,j]/M[i,j])
    end
    return nothing
end

#----------------------------------------------------------------
"""
    normalizerow!(M,i)

Scale row `i` in matrix `M` to make the leading entry 1

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> normalizerow!(M,2)
julia> M
2×3 Array{Rational{Int64},2}:
 1//1  2//1  3//1
 1//1  5//4  3//2
```
"""
function normalizerow!{T<:RealOrSym}(M::Array{T,2},i::Integer)::Void
    ϵ = T <: RatOrInt ? 0 : eps(norm(M,Inf)) 
    j = findfirst(abs(M[i,:]) .> ϵ)
    rowscale!(M,i, T<:RatOrIntOrSym ? 1//M[i,j] : 1/M[i,j])
    return nothing
end


end # module
