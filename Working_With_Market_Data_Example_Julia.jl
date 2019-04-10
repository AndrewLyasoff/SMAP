#=
MIT license (c) 2019 by Andrew Lyasoff

Julia 1.1.0.

The code illustrates some basic operations with large chunks of market data
(extracted from public sources), such as
creating histograms from the returns, and other similar operations. 
=#

using Dates, JuliaDB

#=
The historical quotes are downloaded in .csv format from nasdaq.com and are placed
in the directory "HistoricalQuotes." Before those files can be used here they
must be modified (slightly) in Excel, LibreOffice Calc, or a text editor (e.g., Emacs):
there should be no empty line, or a line that contains a time-stamp (all lines except
the first one must be identically formatted -- this may involve removing the second
line in the spreadhseet). The data for all stocks can be stored in a single variable
(note that HistoricalQuotes is the name of the sub-directory that contains all .csv files).
=#

stocksdata = loadndsparse("HistoricalQuotes"; filenamecol = :ticker, indexcols = [:ticker, :date])

#=
The same database can now be saved in a special binary format that makes reloading
at a later time very fast.
=#
save(stocksdata, "stocksdata.jdb")
@time reloaded_stocksdata = load("stocksdata.jdb")

# Test the saved database for consistency.
stocksdata == reloaded_stocksdata

# Look up the data associated only with Apple (symbol AAPL):
stocksdata["AAPL",:]
length(stocksdata["AAPL",:])

# Similarly, look up the data linked to Microsoft.
stocksdata["MSFT",:]
length(stocksdata["MSFT",:])

# Look up the price of Google on a specific date in the past.
stocksdata["GOOGL", Date(2009,9,9)]
stocksdata["GOOGL", Date(2009,9,9)].close

# Extract the closing prices of Amazon *only*.
selectvalues(stocksdata,:close)["AMZN",:]
keytype(stocksdata)

# List all dates for which the closing prices for Apple are available.
AAPLdate=columns(stocksdata["AAPL",:])[1]

#Extract the closing prices of Apple on those days.
AAPLclose=columns(stocksdata["AAPL",:])[2]

stocksdata["AAPL", Date(2008,8,11)].close

# Produce some plots.
using Plots
pyplot()
plot(AAPLdate, AAPLclose,label="AAPL closing price")


# Calculate and plot the daily returns for AAPL.
begin
    llc=length(AAPLclose);
    AAPLreturns=(AAPLclose[2:llc]-AAPLclose[1:llc-1])./AAPLclose[1:llc-1];
    scatter(AAPLdate[2:llc],AAPLreturns,label="AAPL daily returns",markersize=2)
end


#=
Now build the histogram from the returns using a custom made 'histogram' function
(called 'hstgram' to avoid the confucion with the standard 'histogram').
It takes as an input a single 1-dimensional array of data. The number of bins in
the histogram is determined automatically by using the Diaconis-Friedman rule.
The function returns two arrays: the mid-points of the bins and the (unnormalized)
heights of the bars.
=#

using StatsBase

function hstgram(data_sample::Array{Float64,1})
    data_sorted=sort(data_sample)
    first=data_sorted[1]
    last=data_sorted[end]
    nmb=length(data_sorted)
    IQR=percentile(data_sorted,75)-percentile(data_sorted,25)
    bin_size_loc = 2*IQR*(nmb^(-1.0/3))
    num_bins=Int(floor((last-first)/bin_size_loc))
    bin_size=(last-first)/(num_bins)
    bin_end_points=[first+(i-1)*bin_size for i=1:(num_bins+1)]
    ahist_val=[length(data_sorted[data_sorted .< u]) for u in bin_end_points]
    hist_val=[ahist_val[i+1]-ahist_val[i] for i=1:num_bins]
    mid_bins=[first-bin_size/2+i*bin_size for i=1:num_bins]
    return mid_bins, hist_val
end


# Normalize the bars so that the area of the histogram equals 1.
begin
    U,V=hstgram(AAPLreturns);
    VV=V/(sum(V)*(U[2]-U[1]));
end

begin
    plot(U,VV,line=(:sticks,0.75),label="")
    xlabel!("returns")
    ylabel!("normalized frequency")
    title!("Histogram from the 10y daily returns from AAPL.")
end

# Test that the area of the histogram is indeed 1.
sum(VV)*(U[2]-U[1])


# Another variation of the same histogram.
begin
    plot(U.+(U[2]-U[1])/2,VV,label="",line=(:steppre,1),linewidth=0.05)
    xlabel!("samples")
    ylabel!("normalized frequency")
end

# Normalize the bars to give the probabilities for hitting the bins.
begin
    VVV=V/sum(V);
    plot(U.+(U[2]-U[1])/2,VVV,label="",line=(:steppre,1),linewidth=0.05)
    xlabel!("samples")
    ylabel!("probability")
end

# Test that the probabilities do sum to 1.
sum(VVV)
