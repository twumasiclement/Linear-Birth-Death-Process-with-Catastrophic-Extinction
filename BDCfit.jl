module BDCfit

# for calculating the log likelihood for a BDC model
# uses lambda for birthrate, mu for deathrate, rho for catastrophe rate
# (all per person) and t for time

using PolynomialRoots

export logL
export ProbBDC

function BDCconsts(lambda, mu, rho, t)
    # Constants used in calculating distribution of BDC process at time t
    rts = sort(real(roots([mu, -(lambda+mu+rho), lambda])))
    v0 = rts[1]
    v1 = rts[2]
    sigma = exp(-lambda*(v1 - v0)*t)
    k1 = v0*v1*(1 - sigma)/(v1 - sigma*v0)
    k2 = (v1*sigma - v0)/(v1 - sigma*v0)
    k3 = (1 - sigma)/(v1 - sigma*v0)
    return [k1, k2, k3]
end

function gamma_n_j(nmax)
    # calculates gamma^n_j for n = 1, ..., nmax and j = 1, ..., n
    # used by ProbBDC
    gnj = zeros(BigInt, nmax, nmax)
    gnj[1,1] = 1
    if nmax > 1
        for n = 2:nmax
            gnj[n,1] = n*gnj[n-1,1]
        end
        for j = 2:nmax
            for n = j:nmax
                gnj[n,j] = gnj[n-1,j-1] + (n+j-1)*gnj[n-1,j]
            end
        end
    end
    return gnj
end

function delta_m_j(mmax, k1, k2, k3)
    # calculates delta^m_j for n = 1, ..., mmax and j = 1, ..., n
    # used by ProbBDC; k1, k2, k3 will be output from BDCconsts
    k = (k2 + k1*k3)/k1/k3
    dmj = zeros(BigFloat, mmax, mmax)
    dmj[1,1] = k
    if mmax == 1
        return dmj
    else
        for m = 2:mmax
            dmj[m,1] = k*m
            for j = 2:m
                dmj[m,j] = k*(m - j + 1)*dmj[m,j-1]
            end
        end
        return dmj
    end
end

function ProbBDC(lambda::Float64, mu::Float64, rho::Float64, t::Float64, mmax::Float64, nmax::Float64)
    return ProbBDC(lambda::Float64, mu::Float64, rho::Float64, t::Float64, convert(Int64, mmax), convert(Int64, nmax))
end

function ProbBDC(lambda::Float64, mu::Float64, rho::Float64, t::Int64, mmax::Int64, nmax::Int64)
    return ProbBDC(lambda::Float64, mu::Float64, rho::Float64, convert(Float64, t), mmax, nmax)
end

function ProbBDC(lambda::Float64, mu::Float64, rho::Float64, t::Float64, mmax::Int64, nmax::Int64)
    # P(X_t=n |  X_0=m) for -1 <= m <= mmax and -1 <= n <= nmax
    # where -1 indicates extinction by catastrophe
    cc = BDCconsts(lambda, mu, rho, t)
    k1 = cc[1]
    k2 = cc[2]
    k3 = cc[3]
    k4 = (k1 + k2)/(1 - k3)
    P = zeros(Float64, mmax+2, nmax+2)
    P[1,1] = 1
    P[2,2] = 1
    k1_powers = zeros(BigFloat, mmax)
    k3_powers = zeros(BigFloat, nmax)
    k4_powers = zeros(BigFloat, mmax)
    facts = zeros(BigFloat, nmax)
    k1_powers[1] = k1
    k4_powers[1] = k4
    P[3,1] = Float64(1 - k4)
    P[3,2] = Float64(k1)
    for m = 2:mmax
        k1_powers[m] = k1*k1_powers[m-1]
        k4_powers[m] = k4*k4_powers[m-1]
        P[m+2,1] = Float64(1 - k4_powers[m])
        P[m+2,2] = Float64(k1_powers[m])
    end
    k3_powers[1] = k3
    facts[1] = 1
    for n = 2:nmax
        k3_powers[n] = k3*k3_powers[n-1]
        facts[n] = n*facts[n-1]
    end
    gnj = gamma_n_j(nmax)
    dmj = delta_m_j(mmax, k1, k2, k3)
    for m = 1:mmax
        for n = 1:nmax
            x = BigFloat(0)
            for j = 1:(min(m,n))
                x = x + gnj[n,j]*dmj[m,j]
            end
            P[m+2,n+2] = Float64(x*k1_powers[m]*k3_powers[n]/facts[n])
        end
    end
    return P
end

function logL(lambda::Float64, mu::Float64, rho::Float64, x::Array{Float64, 2})
    return logL(lambda, mu, rho, convert(Array{Int64, 2}, x))
end

function logL(lambda::Float64, mu::Float64, rho::Float64, x::Array{Int64, 2})
    # calculate the log likelihood for params lambda, mu, rho and data x
    # each row of x are population at times 1, 3, 5, 7, 9, 11, 13, 15, 17
    # assume population at time 0 is 2; state -1 indicates catastrophe
    mmax1 = 2
    nmax1 = Int64(max(maximum(x[:,1]), 2))
    P1 = ProbBDC(lambda, mu, rho, 1, mmax1, nmax1)
    mmax2 = Int64(max(maximum(x[:,1:8]), 2))
    nmax2 = Int64(max(maximum(x), 2))
    P2 = ProbBDC(lambda, mu, rho, 2, mmax2, nmax2)
    el = 0
    for i = 1:size(x, 1) # logL for observation i
        # time 0 to time 1 transition
        el = el + log(P1[4, Int64(x[i,1]+2)])
        for j = 1:8
            # time 2j-1 to time 2j+1 transition
            el = el + log(P2[Int64(x[i,j]+2), Int64(x[i,j+1]+2)])
        end
    end
    return el
end

end #module


# using DelimitedFiles
# output_file = "P.csv"
# P = ProbBDC(.3, .2, .001, 1, 1000, 1000)
# writedlm(output_file, P)
# for m = 1:1000
#     print(sum(P[m,:]), "\n")
# end

# BDCdir = "/home/owen/Documents/BDCsim"
# push!(LOAD_PATH, BDCdir)
# using BDCfit
