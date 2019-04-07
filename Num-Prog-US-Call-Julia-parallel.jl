## MIT license (C) 2019 by Andrew Lyasoff
##
## Julia 1.1.0.
## This code implements on multiple CPUs the numerical program for computing the exercise boundary
## and the price of an American-style call option described in Sec. 18.4 of
## "Stochastic Methods in Asset Pricing" (pp. 514-516).
## The same program is implemented in Python in Appendix B.3 of SMAP.
##



begin
    using Distributed
    pnmb=2 #set the number of processes here
    addprocs(pnmb)
    ## The number if grid points in the time domain is set below.
    ## 5000 is excessively large and is meant to test Julia's capabilities.
    total_grid_no=5000 #this number must be divisible by pnmb
    include("Num-Prog-US-Call-Julia-parallel-functions.jl")
    nodes,weights=gausslegendre( 200 );
end

#set the parameters in the model
begin
    KK=40.0; # strike price
    sigma=0.3;
    delta=0.07; # dividend rate
    rr=0.02; # interest rate
    TT=0.5; # time to maturity in years
    loc_no_grd_pnts=Int64(total_grid_no/pnmb); # number of grid points to give to each CPU
    DLT=TT/(total_grid_no); # distance between two grid points in the time domain
    ABSC=range(0.0,length=total_grid_no+1,stop=TT) # the entire grid on the time domain
    VAL=[max(KK,KK*(rr/delta)) for x in ABSC]; # first guess for the exercise boundary
end

## run the main routine
## The third argument is an upper bound on the total iterations to run.
## The last two arguments control FastGaussQuadrature.
## masterF(prnmb,start_iter,nmb_of_iter,conv_tol,K,σ,δ,r,T,Δ,no_grd_pnts,vls,nds,wghts)
##
@time ABSC,VAL,conv,iterations=masterF(pnmb,0,100,1.0e-5,KK,sigma,delta,rr,TT,DLT,loc_no_grd_pnts,VAL,nodes,weights);
## 11.814633 seconds (9.24 M allocations: 486.551 MiB, 2.13% gc time)

# The call to masterF can be repeated with the most recent VAL
# and the second argument set to the number of already performed iterations.

conv #The achived convergence
#9.63110706209136e-6

iterations
#23

f=CubicSplineInterpolation(ABSC,VAL,extrapolation_bc = Interpolations.Line());

begin
    using Plots
    pyplot()
end

begin
    plotgrid=ABSC[1]:.001:ABSC[end];
    pval=[f(x) for x in plotgrid];
    plot(plotgrid,pval,label="exercise boundary")
    xlabel!("real time")
    ylabel!("underlying spot price")
end

f(0)
#53.58270125587188


## The price of the EU at-the-money call with 0.5 years to expiry.
EUcall(KK,TT,KK,sigma,rr,delta)
#2.8378265020118754

## The early exercise premium for a US at-the-money call with 0.5 years to expiry.

begin
    an=[(0.5*(TT-0.0))*exp(-delta*((1/2)*(TT-0.0)*s+(1/2)*(0.0+TT)-0.0))*(KK*delta/2)*F(1,0.0,(1/2)*(TT-0.0)*s+(1/2)*(0.0+TT),f((1/2)*(TT-0.0)*s+(1/2)*(0.0+TT))/KK,rr,delta,sigma) for s in nodes]
    bn=[(0.5*(TT-0.0))*exp(-rr*((1/2)*(TT-0.0)*s+(1/2)*(0.0+TT)-0.0))*(rr*KK/2)*F(-1,0.0,(1/2)*(TT-0.0)*s+(1/2)*(0.0+TT),f((1/2)*(TT-0.0)*s+(1/2)*(0.0+TT))/KK,rr,delta,sigma)  for s in nodes]
    aaa=weights'*an;
    bbb=weights'*bn;
    EEP=aaa-bbb;
end

EEP
#0.10081943488743587

## The early exercise premium for a US at-the-money call with 0.5 years to expiry.
EUcall(KK,TT,KK,sigma,rr,delta)+EEP
#2.938645936899311
