module IntroLinearAlgebra

# package code goes here

export rowswitch, rowscale, rowadd, rowswitch!, rowscale!, rowadd!, rref

import Base.show

RatOrInt = Union{Rational,Integer} 

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

show{T<:RatOrInt}(io::IO, ::MIME"text/latex", s::Array{T,2}) = write(io, texstring(s))

function rowswitch{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer)
    A = copy(M)
    A[i,:], A[j,:] = A[j,:], A[i,:] 
    return A
end

function rowscale{T<:RatOrInt}(M::Array{T,2},i::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
    A[i,:] *= x
    return A
end

function rowadd{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    A = copy(M)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
       A = convert(Array{Rational{Int64},2},A)
    end 
    A[i,:] += x*A[j,:]
    return A
end

function rowswitch!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer)
    if i == j
        return nothing
    else
        M[i,:], M[j,:] = M[j,:], M[i,:]
        return nothing
    end
end

function rowscale!{T<:RatOrInt}(M::Array{T,2},i::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
       	     "to use rowscale! with a rational, like M = [1//1 2; 3 4]")
    end
    M[i,:] *= x
    return nothing
end

function rowadd!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer,x::Real)
    if isa(x,Rational) && issubtype(typeof(M).parameters[1],Integer)
        error("Initialize your matrix with at least one rational number " *
       	     "to use rowadd! with a rational, like M = [1//1 2; 3 4]")
    else
        M[i,:] += x*M[j,:]
    end
    return nothing
end


function rref{T<:RatOrInt}(M::Array{T,2};showsteps=false)
    A = convert(Array{Rational{Int64},2},M)
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
        return steps
    else
        return A
    end
end

function zerooutbelow!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer) # zero out entries below the (i,j)th
    for k=i+1:size(M,1)
        rowadd!(M,k,j,-M[k,j]//M[i,j])
    end
    return nothing
end


function zerooutabove!{T<:RatOrInt}(M::Array{T,2},i::Integer,j::Integer) 
    for k=1:i-1
        rowadd!(M,k,j,-M[k,j]//M[i,j])
    end
    return nothing
end


function normalizerow!{T<:RatOrInt}(M::Array{T,2},i::Integer)
    j = findfirst(M[i,:] .!= 0)
    rowscale!(M,i,1//M[i,j])
    return nothing
end


end # module
