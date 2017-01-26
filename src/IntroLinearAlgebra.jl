module IntroLinearAlgebra

# package code goes here

export rowswitch, rowscale, rowadd

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

end # module
