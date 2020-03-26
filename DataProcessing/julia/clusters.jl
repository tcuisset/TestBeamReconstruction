using DataFrames
using CSV
using Query
using Plots
plotly()
#functions
showln(x) = (show(x); println())
log_weight(weight, total) = params.thresholdW0 + log(weight/total)

#declare an immutable structur to hold the parameters
#structures in Julia do not contain methods, only variables
#declaring a variable without a type is equivalent to var::Any
struct Params
    thresholdW0::Float32
end
params = Params(2.9)

#strings are added with '*' instead of '+'!
#'.' is the element-wise operator
basename = "../Clue/data/output/cluster_data_"
filenames = basename .* ("20_5_20_20_0.csv",
                         "20_5_20_20_1.csv")

print(filenames)
for file in filenames
    #convert the input CSV files into DataFrames objects
    df_init = CSV.read(file)

    #LINQ-like query sintax (it feels like SQL)
    #I could have also used the native DataFrames selection tools or
    #the standalone query operators of the Query package
    #http://www.queryverse.org/Query.jl/stable/linqquerycommands/
    #Merging multiple queries into one does not seem to be supported
    df1 = @from i in df_init begin
        @group i by i.clusterId into g
        @select {ClusterId=key(g), NHits=length(g), WeightSum=sum(g.weight)}
        @collect DataFrame
    end

    df2 = @from i in df_init begin
        @join j in df1 on i.clusterId equals j.ClusterId
        @select {i.clusterId, i.weight, xmult=i.x*log_weight(i.weight,j.WeightSum), ymult=i.y*log_weight(i.weight,j.WeightSum)}
        @collect DataFrame
    end

    df3 = @from i in df2 begin
        @group i by i.clusterId into g
        @select {ClusterId=key(g), XLog=sum(g.xmult), YLog=sum(g.ymult), WeightSum=sum(g.weight)}
        @collect DataFrame
    end

    #see Arabella's calculatePositions() code: 
    #https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalCLUEAlgo.cc#L183-L237
    xlog = convert(Array{Float32}, df3[!,:XLog])
    ylog = convert(Array{Float32}, df3[!,:YLog])
    weightsum = convert(Array{Float32}, df3[!,:WeightSum])
    xlog_weighted = xlog ./ weightsum
    ylog_weighted = ylog ./ weightsum

    println("Cluster Positions:")
    println("x: ", xlog_weighted)
    println("y: ", ylog_weighted)

    p = plot(xlog_weighted, ylog_weighted, seriestype=:scatter)
    savefig(p, "test.png")
end #for file in filenames:


#=
import numpy as np
import pandas as pd

dataframes = ( pd.read_csv('../Clue/data/output/cluster_data_20_5_20_20_0.csv'),
               pd.read_csv('../Clue/data/output/cluster_data_20_5_20_20_1.csv') )

def plot_energy():
    for df in dataframes:
        mask_clusters = (df['clusterId'] != -1)
        mask_outliers = (df['clusterId'] == -1)
        clusters = df.loc[mask_clusters]
        outliers = df.loc[mask_outliers]

        en_tot = df.loc[:,'weight'].sum()
        en_clusters = clusters.loc[:,'weight'].sum()
        en_outliers = outliers.loc[:,'weight'].sum()
        en_outliers_fraction = en_outliers / en_tot
        print(en_tot, en_clusters, en_outliers, en_outliers_fraction)

        nhits = clusters.groupby('clusterId').size()
        print(nhits)

def plot_position():
    for df in dataframes:
        pass
    

def main():
    #plot_energy()
    plot_position()

if __name__ == '__main__':
    main()
=#
