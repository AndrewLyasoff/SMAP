## MIT license (c) 2019 by Andrew Lyasoff
##
## note the formula
## $\int_a^bf(t)dt={1\over2}(b-a)\int_{-1}^{+1}f\bigl({1\over2}(b-a)s+{1\over2}(a+b)\bigr)d s$
##(Used by FastGaussQuadrature)
##



begin
    @everywhere using SpecialFunctions
    @everywhere using Interpolations
    @everywhere using Roots
    @everywhere using FastGaussQuadrature

    @everywhere function EUcall(S::Float64,t::Float64,K::Float64,σ::Float64,r::Float64,δ::Float64)
        local v1,v2
        v1=(δ-r+σ^2/2.0)*t+log(K/S)
        v2=(-δ + r + σ^2/2.0)*t - log(K/S)
        v3=σ*sqrt(2.0*t)
        return -K*(exp(-t*r)/2)+K*(exp(-t*r)/2)*erf(v1/v3)+S*(exp(-t*δ)/2)*erf(v2/v3)+S*(exp(-t*δ)/2)
    end

    @everywhere function F(ϵ::Int64,t::Float64,u::Float64,v::Float64,r::Float64,δ::Float64,σ::Float64)
        v1=(r-δ+ϵ*σ^2/2)*(u-t)-log(v)
        v2=σ*sqrt(2*(u-t))
        return 1.0+erf(v1/v2)
    end

    @everywhere function ah(t::Float64,z::Float64,r::Float64,δ::Float64,σ::Float64,f)
        return (exp(-δ*(z-t))*(δ/2)*F(1,t,z,f(z)/f(t),r,δ,σ))
    end

    @everywhere function bh(t::Float64,z::Float64,r::Float64,δ::Float64,σ::Float64,K::Float64,f)
        return (exp(-r*(z-t))*(r*K/2)*F(-1,t,z,f(z)/f(t),r,δ,σ))
    end

    @everywhere function mke_grd(bgn::Float64,step::Float64,size::Int64)
        return bgn:step:(bgn+step+(size-1)*step)
    end


    @everywhere function workerF(loc_range::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},r::Float64,δ::Float64,σ::Float64,T::Float64,K::Float64,nds::Array{Float64,1},wghts::Array{Float64,1},f)
        local loc,LRT,aaa,bbb,an,bn,lnlr
        lnlr=length(loc_range)
        loc=zeros(lnlr)
        for i=lnlr:-1:1
            t=loc_range[i]
            an=[(1/2)*(T-t)*ah(t,(1/2)*(T-t)*s+(1/2)*(t+T),r,δ,σ,f) for s in nds]
            bn=[(1/2)*(T-t)*bh(t,(1/2)*(T-t)*s+(1/2)*(t+T),r,δ,σ,K,f) for s in nds]
            aaa=wghts'*an
            bbb=wghts'*bn;
            LRT=find_zero(x->x-K-EUcall(x,T-t,K,σ,r,δ)-aaa*x+bbb,(K-10,K+20));
            loc[i]=LRT
        end
        return loc
    end

    function cmbn(varg::Array{Array{Float64,1},1})
        if length(varg)>1
            loc=varg[1]
            for i=2:length(varg)
                loc=vcat(loc,varg[i])
            end
            return loc
        else
            return varg
        end
    end

end


function masterF(prnmb::Int64, start_iter::Int64,nmb_of_iter::Int64,conv_tol::Float64,K::Float64,σ::Float64,δ::Float64,r::Float64,T::Float64,Δ::Float64,no_grd_pnts::Int64,vls::Array{Float64,1},nds::Array{Float64,1},wghts::Array{Float64,1})
    local no_iter,conv_check,absc,absc0,val,valPrev,loc,f
    absc=range(0.0,length=pnmb*no_grd_pnts+1,stop=T)
    absc0=range(0.0,length=pnmb*no_grd_pnts,stop=(T-Δ))
    all_grid=[range((T/prnmb)*i,length=no_grd_pnts,stop=((T/prnmb)*(i+1)-Δ)) for i=0:(pnmb-1)]
    val=vls;
    f=CubicSplineInterpolation(absc,val,extrapolation_bc = Interpolations.Line())
    for i=2:(prnmb+1)
        remotecall_fetch(()->K, i)
        remotecall_fetch(()->σ, i)
        remotecall_fetch(()->δ, i)
        remotecall_fetch(()->r, i)
        remotecall_fetch(()->T, i)
        remotecall_fetch(()->all_grid, i)
        remotecall_fetch(()->f, i)
        remotecall_fetch(()->nds, i)
        remotecall_fetch(()->wghts, i)
    end
    no_iter=start_iter;
    conv_check=1000.0
    while no_iter<nmb_of_iter&&conv_check>conv_tol
        valPrev=val;
        no_iter+=1
        lst=[@spawnat i workerF(all_grid[myid()-1],r,δ,σ,T,K,nds,wghts,f) for i=2:(prnmb+1)]
        val0=[fetch(lst[i]) for i=1:prnmb];
        val=cmbn(val0)
        push!(val,max(K,K*(r/δ)))
        f=CubicSplineInterpolation(absc,val,extrapolation_bc = Interpolations.Line())
        for i=2:(prnmb+1)
            remotecall_fetch(()->f, i)
        end
        conv_check=maximum(abs.(valPrev-val))
    end
    return absc,val,conv_check,no_iter
end
