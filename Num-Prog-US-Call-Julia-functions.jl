## MIT license (C) 2019 by Andrew Lyasoff
##
## note the formula
## $\int_a^bf(t)dt={1\over2}(b-a)\int_{-1}^{+1}f\bigl({1\over2}(b-a)s+{1\over2}(a+b)\bigr)d s$
##(Used by FastGaussQuadrature)
##




begin
    using SpecialFunctions
    using Interpolations
    using Roots
    using FastGaussQuadrature


    function EUcall(S::Float64,t::Float64,K::Float64,σ::Float64,r::Float64,δ::Float64)
        local v1,v2
        v1=(δ-r+σ^2/2.0)*t+log(K/S)
        v2=(-δ + r + σ^2/2.0)*t - log(K/S)
        v3=σ*sqrt(2.0*t)
        return -K*(exp(-t*r)/2)+K*(exp(-t*r)/2)*erf(v1/v3)+S*(exp(-t*δ)/2)*erf(v2/v3)+S*(exp(-t*δ)/2)
    end

    function F(ϵ::Int64,t::Float64,u::Float64,v::Float64,r::Float64,δ::Float64,σ::Float64)
        v1=(r-δ+ϵ*σ^2/2)*(u-t)-log(v)
        v2=σ*sqrt(2*(u-t))
        return 1.0+erf(v1/v2)
    end

    function ah(t::Float64,z::Float64,r::Float64,δ::Float64,σ::Float64,f)
        return (exp(-δ*(z-t))*(δ/2)*F(1,t,z,f(z)/f(t),r,δ,σ))
    end

    function bh(t::Float64,z::Float64,r::Float64,δ::Float64,σ::Float64,K::Float64,f)
        return (exp(-r*(z-t))*(r*K/2)*F(-1,t,z,f(z)/f(t),r,δ,σ))
    end

    function make_grid0(step::Float64,size::Int64)
        return 0.0:step:(step+(size-1)*step)
    end
end

function mainF(start_iter::Int64,nmb_of_iter::Int64,conv_tol::Float64,K::Float64,σ::Float64,δ::Float64,r::Float64,T::Float64,Δ::Float64,nmb_grd_pnts::Int64,vls::Array{Float64,1},nds::Array{Float64,1},wghts::Array{Float64,1})
    local no_iter,conv_check,absc,absc0,valPrev,val,loc,f
    absc=range(0.0,length=nmb_grd_pnts+1,stop=T)
    absc0=range(0.0,length=nmb_grd_pnts,stop=(T-Δ));
    val=vls;
    f=CubicSplineInterpolation(absc,val,extrapolation_bc = Interpolations.Line())
    no_iter=start_iter;
    conv_check=100.0
    while no_iter<nmb_of_iter&&conv_check>conv_tol
        no_iter+=1
        loc=[max(K,K*(r/δ))]
        for ttt=Iterators.reverse(absc0)
            an=[(1/2)*(T-ttt)*ah(ttt,(1/2)*(T-ttt)*s+(1/2)*(ttt+T),r,δ,σ,f) for s in nodes]
            bn=[(1/2)*(T-ttt)*bh(ttt,(1/2)*(T-ttt)*s+(1/2)*(ttt+T),r,δ,σ,K,f) for s in nodes]
            aaa=weights'*an
            bbb=weights'*bn;
            LRT=find_zero(x->x-K-EUcall(x,T-ttt,K,σ,r,δ)-aaa*x+bbb,(K-10,K+20));
            pushfirst!(loc,LRT)
        end
        valPrev=val
        val=loc
        f=CubicSplineInterpolation(absc,val,extrapolation_bc = Interpolations.Line())
        conv_check=maximum(abs.(valPrev-val))
    end
    return absc,val,conv_check,no_iter
end
