module GaussLegendreRKRule

using LinearAlgebra
include("LegendrePoly.jl")

export buildRule

struct ButcherTable
    num_stages::Int
    c
    A
    b
    btilde
end

function buildRule(n::Int; outputType = Float64)

    # determine zeros of shifted Legendre polynomial,on [0,1]
    c = (LegendrePoly.getZeros(n) .+ big"1")./big"2"

    # build Vandermonde matrix
    V = hcat( ones(BigFloat, n, 1), [ c[i]^(j-1) for i = 1:n, j = 2:n ] )

    # build $e_H$
    eH = [ big"1"/i for i in 1:n ]

    # build $\tilde{e}_H$, reduced order s-1
    etH = deepcopy(eH)
    etH[end] = 0;

    # Build A
    C = [ c[i]^j/j for i in 1:n, j in 1:n ]
    A = C/V

    # Build full order b
    b = (eH'/V)'

    # build reduced order bt
    bt = (etH'/V)'

    #return ( c, A, b, bt )
    if outputType == Float64
        return ButcherTable( n,
            Float64.(c),
            Float64.(A),
            Float64.(b),
            Float64.(bt) )
    elseif outputType == BigFloat
        return ButcherTable( n, c, A, b, bt )
    else
        println("Error: requested type has not been implemented")
    end

end


end
