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
       normalizerow!

import Base.show

if Pkg.installed("SymPy") != nothing
    import SymPy
    sympyexists = true
else 
    sympyexists = false
end

RatOrInt = sympyexists ? Union{Rational,Integer,SymPy.Sym} : Union{Rational,Integer}

"""
    textstring(M)

Return a string for processing by LaTeX to pretty print
the matrix `M`
"""
function texstring{T<:RatOrInt}(M::Array{T,2})
    function prettystring(x::RatOrInt)
        if isa(x,Integer)
            return string(x)
        elseif x.den == 1
            return string(x.num)
        else
            sgn = signbit(x.num) ? "-" : ""
            return sgn*"\\frac{"*string(abs(x.num))*"}{"*string(x.den)*"}"
        end
    end
    s = "\$\\left[\\begin{array}{" * repeat("c",size(M,2)) * "}"
    for i=1:size(M,1)
        for j=1:size(M,2)
            s *= prettystring(M[i,j]) 
            if j < size(M,2)
                s *= " & "
            end
        end
        s *= " \\\\ "
    end
    s *= "\\end{array}\\right]\$"
    return s
end

if !sympyexists
    show{T<:RatOrInt}(io::IO, ::MIME"text/latex", s::Array{T,2}) = write(io, texstring(s))
end

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
function rowswitch{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer)
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
function rowscale{T<:RatOrInt}(M::Array{T,2},i::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
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
function rowadd{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
    A[i,:] += x*A[j,:]
    return A
end

"""
    rowswitch!(M,i,j)

Version of `rowswitch` that modifies the argument rather than 
returning a new matrix
"""
function rowswitch!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer)
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
function rowscale!{T<:RatOrInt}(M::Array{T,2},i::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
                    "to use rowscale! with a rational, like M = [1//1 2; 3 4]")
    end
    M[i,:] *= x
    return nothing
end

"""
    rowadd!(M,i,j)

Version of `rowadd` that modifies the argument rather than 
returning a new matrix
"""
function rowadd!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
                    "to use rowadd! with a rational, like M = [1//1 2; 3 4]")
    else
        M[i,:] += x*M[j,:]
    end
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
function rref{T<:RatOrInt}(M::Array{T,2};showsteps=false)
    A = copy(convert(Array{Rational{Int64},2},M))
    steps = typeof(A)[]
    current_row = 1
    for j=1:size(A,2)
        if all(A[:,j] .== 0)
            continue
        else
            i = current_row-1 + findfirst(A[current_row:end,j] .!= 0)
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
        if all(A[i,:] .== 0)
            continue
        else
            j = findfirst(A[i,:] .!= 0) 
            zerooutabove!(A,i,j)
            if showsteps push!(steps,copy(A)) end
        end
    end
    if showsteps
        uniquesteps = Array{Rational{Int64},2}[]
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
function zerooutbelow!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer) # zero out entries below the (i,j)th
    for k=i+1:size(M,1)
        rowadd!(M,k,i,-M[k,j]//M[i,j])
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
function zerooutabove!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer)::Void
    for k=1:i-1
        rowadd!(M,k,i,-M[k,j]//M[i,j])
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
function normalizerow!{T<:RatOrInt}(M::Array{T,2},i::Integer)::Void
    j = findfirst(M[i,:] .!= 0)
    rowscale!(M,i,1//M[i,j])
    return nothing
end


end # module
