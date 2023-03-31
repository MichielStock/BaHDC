#=
Created on 31/03/2023 13:47:33
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Analysing the odor dataset.
=#


using DrWatson
@quickactivate "BaHDC"
using DataFrames, CSV, LinearAlgebra

# loading the data
# ----------------


odors = CSV.read(datadir("odor/TrainSet.txt"), DataFrame)

# features
dragon = CSV.read(datadir("odor/Dragon"), DataFrame)
fingerprints = CSV.read(datadir("odor/Fingerprints"), DataFrame)
descriptors = CSV.read(datadir("odor/molecular_descriptors_data.txt"), DataFrame)

smiles = CSV.read(datadir("odor/Smiles"), DataFrame)


# analysis
# --------

# make some plots about the distribution of the data


# you can transform features in a HDV embedding using a random mapping
# CAUTION: you might need to normalize the columns....

features = randn(190, 100)  # feature matrix

W = randn(100, 10_000)  # random mapping to 10_000 dims

X = sign.(features * W) .|> Int


# embed the odor profile
# perform SVD