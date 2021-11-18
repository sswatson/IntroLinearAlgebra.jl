__precompile__(true)

module IntroLinearAlgebra

using Pkg, LinearAlgebra, Requires

export rowswitch,
       rowscale,
       rowadd,
       rref,
       eigenvalues,
       eigenspaces,
       charpoly,
       zerooutabove!,
       zerooutbelow!,
       rownormalize!,
       setdigits,
       setalignment,
       transformation_movie,
       svdviz

import Base: show, eps, inv

#-------------------------------------------------------
# MATRIX TOOLS
#-------------------------------------------------------

DIGITS = 4
ALIGNMENT = true

RatOrInt = Union{Rational,Integer}

"""
    setdigits!(n::Integer)

Set the number of digits (after the decimal point) to 
display (for matrices with floating-point entries). 
"""
function setdigits!(n::Integer)
    global DIGITS
    DIGITS = n
end

"""
    setalignment!(b::Bool)

Determine whether to align displayed matrix entries
disregarding the sign (true) or including the sign (false)
"""
function setalignment!(b::Bool)
    global ALIGNMENT
    ALIGNMENT = b
end

"""
    texstring(M)

Return a string for processing by LaTeX to pretty print
the matrix `M`
"""
function texstring(M::AbstractArray{T,2}) where T 
    ϵ = T<:RatOrInt ? 0 : eps(norm(M,Inf))
    function neg(x::Real) # test whether a number is
                          # genuinely negative (-0.0 isn't)
        return x < -ϵ
    end
    prettystring(x::Bool;padding=false) = "\\mathrm{"*(x ? "true" : "false")*"}"
    function prettystring(x::RatOrInt;padding=false)
        if isa(x,Integer)
            if padding
                return neg(x) ? string(x) : "\\hphantom{-}"*string(x)
            else
                return string(x)
            end
        elseif x.den == 1
            return prettystring(x.num,padding=padding)
        else
            sgn = neg(x.num) ? "-" : (padding ? "\\hphantom{-}" : "")
            return sgn*"\\frac{"*string(abs(x.num))*"}{"*string(x.den)*"}"
        end
    end
    function prettystring(x::Real;padding=false)
        global DIGITS
        return (neg(x) ? "-" : (padding ? "\\hphantom{-}" : "")) *
            string(round(abs(x);digits=DIGITS))
    end
    function prettystring(z::Complex;padding=false)
        if ~_isnonzero(z.im,[one(typeof(z.im))])
            return prettystring(z.re)
        elseif ~_isnonzero(z.re,[one(typeof(z.re))])
            return (abs(z.im) == 1 ? "" : prettystring(z.im))*"i"
        else
            return prettystring(z.re)*(neg(z.im) ? "-" : "+")*prettystring(abs(z.im)*im)
        end
    end

    s = "\$\\left[\\begin{array}{" * repeat("c",size(M,2)) * "}"
    global ALIGNMENT
    columnpadding = T <: Real ? any(M .< -ϵ, dims=1) : repeat([false],size(M,2))

    for i=1:size(M,1)
        for j=1:size(M,2)
            thiscolumnpadding = T <: Real ? any(M[:,j] .< -ϵ) : false
            s *= prettystring(M[i,j],padding=ALIGNMENT&&columnpadding[j])
            if j < size(M,2)
                s *= " & "
            end
        end
        s *= " \\\\ "
    end
    s *= "\\end{array}\\right]\$"
    return s
end

texstring(M::AbstractArray{T,1}) where T<:Real = texstring(reshape(M,(length(M),1)))
texstring(M::Adjoint{T, Vector{T}}) where T<:Real = texstring(convert(Array{T,2},M))

show(io::IO,
     ::MIME"text/latex",
     s::AbstractArray{T,2}) where T <: Real = write(io, texstring(s))

show(io::IO,
     ::MIME"text/latex",
     s::AbstractArray{T,1}) where T <: Real = write(io, texstring(s))

show(io::IO,
     ::MIME"text/latex",
     s::AbstractArray{T,2}) where T <: Complex = write(io, texstring(s))

show(io::IO,
     ::MIME"text/latex",
     s::AbstractArray{T,1}) where T <: Complex = write(io, texstring(s))

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
function rowswitch(M::AbstractArray{T,2},i::Integer,j::Integer) where T 
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
function rowscale(M::AbstractArray{T,2},i::Integer,x) where T 
    A = copy(convert(Array{promote_type(T,typeof(x)),2},M))
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
function rowadd(M::AbstractArray{T,2},
                i::Integer,
                j::Integer,
                x) where T 
    A = copy(convert(Array{promote_type(T,typeof(x)),2},M))
    A[i,:] += x*A[j,:]
    return A
end

"""
    rowswitch!(M,i,j)

Version of `rowswitch` that modifies the argument rather than 
returning a new matrix
"""
function rowswitch!(M::AbstractArray{T,2},i::Integer,j::Integer) where T 
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
function rowscale!(M::AbstractArray{T,2},i::Integer,x) where T 
    M[i,:] *= x
    return nothing
end

"""
    rowadd!(M,i,j)

Version of `rowadd` that modifies the argument rather than 
returning a new matrix
"""
function rowadd!(M::AbstractArray{T,2},i::Integer,j::Integer,x) where T 
    M[i,:] += x*M[j,:]
    return nothing
end

"""
    _isnonzero(x,M)

Check whether x is nonzero. Symbols count as nonzero, 
and floats count as zero whenever their absolute value
is less than roundoff error for floats the size of the 
norm of M. 
"""
function _isnonzero(x,M::AbstractArray{T}) where T 
    if isa(x,AbstractFloat) || (isa(x,Complex) && isa(x.re,AbstractFloat))
        ϵ = T <: AbstractFloat ? eps(norm(M,Inf)) : eps(typeof(x))
        return abs(x) > ϵ
    else
        return x != 0
    end
end
eps(::Type{Complex{T}}) where T = 3*eps(T)
eps(z::Complex{T}) where T <: AbstractFloat = 3*(eps(z.re) + eps(z.im))

simplify(x) = isa(x,Real) ? x :
    (string(typeof(x)) == "SymPy.Sym" ? symplify(x) : x)

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
function rref(M::AbstractArray{T,2};showsteps=false) where T 
    isnonzero(x) = _isnonzero(x,M)
    S = T <: Rational ? T.parameters[1] : T
    if T <: RatOrInt
        newtype = Rational{S}
    elseif T <: Complex && T.parameters[1] <: Integer
        newtype = Complex{Rational{T.parameters[1]}}
    elseif T <: Complex && T.parameters[1] <: Rational
        newtype = Complex{T.parameters[1]}
    else
        newtype = promote_type(typeof(inv(one(T))),T)
    end
    A = copy(convert(Array{newtype,2},M))
    steps = typeof(A)[]
    current_row = 1
    for j=1:size(A,2)
        if ~any(isnonzero.(A[current_row:end,j]))
            if issubtype(T,AbstractFloat)
                A[current_row:end,j] = 0
            end
        else
            i = current_row-1 + findfirst(isnonzero.(A[current_row:end,j]))
            rowswitch!(A,current_row,i)
            if showsteps push!(steps,copy(A)) end
            zerooutbelow!(A,current_row,j)
            if showsteps push!(steps,copy(A)) end
            rownormalize!(A,current_row)
            if showsteps push!(steps,copy(A)) end
            current_row += 1
            if current_row > size(A,1)
                break
            end
        end
        map!(simplify,A,A)
    end
    for i=size(A,1):-1:1
        if ~any(isnonzero.(A[i,:]))
            continue
        else
            j = findfirst(isnonzero.(A[i,:]))
            zerooutabove!(A,i,j)
            if showsteps push!(steps,copy(A)) end
        end
    end
    if showsteps
        uniquesteps = typeof(A)[]
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
function zerooutbelow!(M::Array{T,2},i::Integer,j::Integer)::Nothing where T
    for k=i+1:size(M,1)
        rowadd!(M,k,i,T<:RatOrInt ? -M[k,j]//M[i,j] : -M[k,j]/M[i,j])
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
function zerooutabove!(M::Array{T,2},i::Integer,j::Integer)::Nothing where T 
    for k=1:i-1
        rowadd!(M,k,i,T<:RatOrInt ? -M[k,j]//M[i,j] : -M[k,j]/M[i,j])
    end
    return nothing
end

#----------------------------------------------------------------
"""
    rownormalize!(M,i)

Scale row `i` in matrix `M` to make the leading entry 1

```julia 
julia> M = [1//1 2 3; 4 5 6]
julia> rownormalize!(M,2)
julia> M
2×3 Array{Rational{Int64},2}:
 1//1  2//1  3//1
 1//1  5//4  3//2
```
"""
function rownormalize!(M::Array{T,2},i::Integer)::Nothing where T
    isnonzero(x) = _isnonzero(x,M)
    j = findfirst(isnonzero.(M[i,:]))
    rowscale!(M,i, T<:RatOrInt ? 1//M[i,j] : 1/M[i,j])
    return nothing
end

function inv(A::Union{Array{Rational{T},2},Array{T,2}}) where T <: Integer
    if size(A,1) != size(A,2)
        error("Matrix not square")
    end
    n = size(A,1)
    E = rref([A I])
    C,D = E[:,1:n], E[:,n+1:end]
    if C != I
        error("Not invertible")
    else
        return D
    end
end


#----------------------------------------------------------------------------
# GRAPHICS TOOLS
#----------------------------------------------------------------------------


function __init__()
    
    Requires.@require SymPy="24249f21-da20-56a4-8eb1-6a02cf4ae2e6" begin

        symplify(x) = SymPy.simplify(x)

        texstring(A::Array{SymPy.Sym}) = SymPy.sympy.latex(A)

        struct EigenSpace
            λ
            mult
            eigenmatrix
        end
        
        function texstring(E::EigenSpace)
            s = "\\mathrm{eigenspace}\\left("
            s *= "λ = "
            s *= SymPy.latex(E.λ,itex=true)*", "
            s *= "\\mathrm{mult} = "*string(E.mult)*", "
            s *= "\\mathrm{Col}\\left("
            io = IOBuffer()
            show(io,(MIME"text/latex").instance,E.eigenmatrix)
            s *= String(take!(io))
            s *= "\\right)\\right)"
        end
        
        function show(io::IO,M::MIME"text/latex",E::EigenSpace)
            write(io,"\$"*texstring(E)*"\$")
        end

        function show(io::IO,E::EigenSpace)
            s = "eigenspace(\n"
            s *= "λ = "
            _io = IOBuffer()
            show(_io,(MIME"text/plain").instance,E.λ)
            s *= String(take!(_io))
            s *= ", mult = "*string(E.mult)*", "
            show(_io,E.eigenmatrix)
            s *= String(take!(_io))
            s *= "\n)"
            print(io,s)
        end

        function show(io::IO,::MIME"text/plain",A::Array{EigenSpace,1})
            for (k,E) in enumerate(A)
                show(io,E)
                if k < length(A)
                    print(io,"\n\n")
                end
            end
        end

        function show(io::IO,M::MIME"text/latex",A::Array{EigenSpace,1})
            print(io,"\$\$\\begin{bmatrix}")
            for E in A
                write(io,texstring(E))
                print(io,"\\\\")
            end
            print(io,"\\end{bmatrix}\$\$")
        end

        """
        eigenvalues(M)

        Uses SymPy to find the exact eigenvalues of M

        ```julia 
        julia> M = [1 2; 4 5]
        julia> eigenvalues(M)
        2-element Array{SymPy.Sym,1}
        ⎡3 + 2⋅√3 ⎤
        ⎢         ⎥
        ⎣-2⋅√3 + 3⎦
        ```
        """
        function eigenvalues(A::Array{T,2}) where T 
            λ = SymPy.symbols("lambda")
            return SymPy.solve(charpoly(A), λ)
        end

        """
        eigenvalues(M)

        Use SymPy to find the exact eigenvalues of M

        ```julia 
        julia> M = [1 2; 4 5]
        julia> eigenspaces(M)
        eigenspace(
        λ = 3 + 2 \\sqrt{3}, mult = 1, 
        ⎡  1   √3⎤
        ⎢- ─ + ──⎥
        ⎢  2   2 ⎥
        ⎢        ⎥
        ⎣   1    ⎦
        )

        eigenspace(
        λ = - 2 \\sqrt{3} + 3, mult = 1, 
        ⎡  √3   1⎤
        ⎢- ── - ─⎥
        ⎢  2    2⎥
        ⎢        ⎥
        ⎣   1    ⎦
        )
        ```
        """
        function eigenspaces(A::Array{T,2}) where T 
            B = map(x->convert(SymPy.Sym,x),A)
            return [EigenSpace(a,b,map(SymPy.simplify,hcat(c...))) for
                    (a,b,c) in B.eigenvects()]
        end

        """
        charpoly(M)

        Uses SymPy to find the characteristic polynomial of M
    
        ```julia 
        julia> M = [1 2; 4 5]
        julia> charpoly(M)
        ```
        """
        function charpoly(A::Array{T,2}) where T 
            B = map(x->convert(SymPy.Sym,x), A)
            return B.charpoly()
        end
        
    end

    Requires.@require AsyPlots="77e5a97a-5ef9-58df-9d21-21957d92d960" begin 
        
        """
        transformation_movie(A,n;frames=20)

        Return an array of Graphics2D arrays which give a dynamic
        picture of how A transforms integer grid lines in [-n,n]^2

        ```julia
        julia> A = [0 2; -1 0]
        julia> movie = transformation_movie(A,10,20);
        julia> using Interact
        julia> @manipulate for i=1:length(movie)
                   movie[i]
               end
        ```
        """
        transformation_movie(A::Array{T,2};kwargs...) where T =
            transformation_movie((x,y)->A*[x;y];kwargs...)

        function transformation_movie(f::Function;
                                      gridlines::Integer=8,
                                      frames::Integer=20)

            n = gridlines

            Path = AsyPlots.Path
            Point = AsyPlots.Point
            NamedColor = AsyPlots.NamedColor

            mypath(A::Array{<:Real,2};kwargs...) = Path([tuple(A[k,:]...) for k=1:size(A,1)];kwargs...)
            mypath(A::Array{T,1};kwargs...) where T = Path([tuple(v...) for v in A];kwargs...)
            sc(r::Real,s::String) = r*NamedColor(s)

            m = 2.05*max(norm(f(n,n),Inf),norm(f(-n,n),Inf))
            box = mypath([-m -m; m -m; m m; -m m; -m -m];color="white")
            return [AsyPlots.Plot(vcat([[Path([0 0; m 0]),Path([0 0; m 0])];
                                        [mypath([(1-t)*[j,k] + t*f(j,k) for j=-n:0.1:n];
                                                linewidth=1.5,color=sc(0.9,"DarkRed")) for k=-n:n];
                                        [mypath([(1-t)*[j,k] + t*f(j,k) for k=-n:0.1:n];
                                                linewidth=1.5,color=sc(0.9,"MidnightBlue")) for j=-n:n];
                                        [Point(0,0;linewidth=1e-2)];
                                        [box]]...))
                    for t=range(0,stop=1,length=frames)]
        end

        function svdviz(A::Array{T,2},
                        B::Array{T,2}=A;
                        n=5,
                        frames=20,
                        gridlinewidth=1,
                        graylinewidth=1,
                        pointsize=1,
                        kwargs...) where T <: Real

            Path = AsyPlots.Path
            Point = AsyPlots.Point
            Arrow = AsyPlots.Arrow
            Plot = AsyPlots.Plot
            
            U, Σ, V = svd(A)
            m = 1.05*max(norm(A*V*[n,n],Inf),norm(A*V*[-n,n],Inf),norm(A*V*[-n,-n],Inf))
            mybox = Path([-m -m; m -m; m m; -m m; -m -m];color="white")
            loc(j,k,t) = tuple((j*B*V[:,1] + k*B*V[:,2])*t + (j*V[:,1] + k*V[:,2])*(1-t)...)
            pts(t) = vcat([Point(loc(j,k,t);linewidth=pointsize) for j=-n:n,k=-n:n]...)
            bluelines(t) = [Path([loc(j,k,t) for j=-n:n];linewidth=gridlinewidth,color="MidnightBlue") for k=-n:n]
            redlines(t) = [Path([loc(j,k,t) for k=-n:n];linewidth=gridlinewidth,color="DarkRed") for j=-n:n]
            graylines = [[Path([(j,k) for j=-n:n];color="Gray",linewidth=graylinewidth,opacity=0.3) for k=-n:n];
                         [Path([(j,k) for k=-n:n];color="Gray",linewidth=graylinewidth,opacity=0.3) for j=-n:n]]
            axes = [Path([0 0; 0 m]),Path([0 0; m 0])]
            basisvectors(t) = [Path([(0,0),tuple([1,0]*(1-t) + A*[1,0]*(t)...)];arrow=Arrow(),linewidth=1.5,color="SeaGreen"),
                               Path([(0,0),tuple([0,1]*(1-t) + A*[0,1]*(t)...)];arrow=Arrow(),linewidth=1.5,color="DarkRed")]
            [Plot([graylines;
                   basisvectors(t);
                   bluelines(t);
                   redlines(t);
                   pts(t);
                   [mybox]];kwargs...) for t = range(0,stop=1,length=frames)]
        end

    end

end

end # module
