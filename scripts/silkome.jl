#=
Created on 30/03/2023 12:25:56
Last update: 31/03/2023

@author: Michiel Stock
michielfmstock@gmail.com

Loader for the rough silk data loader.

Data from 

K. Arakawa et al., “1000 spider silkomes: Linking sequences to silk physical properties,” Science Advances, vol. 8, no. 41, p. eabo6043, Oct. 2022, doi: 10.1126/sciadv.abo6043.

=#



using DrWatson
@quickactivate "BaHDC"
using DataFrames, FASTX, CSV, LinearAlgebra

# loading the data
# ----------------

# load organsims 🕷

organisms = CSV.read(datadir("silkome/organisms.csv"), DataFrame)

org2id = Dict(zip(zip(organisms[:,:family], organisms[:,:genus], organisms[:,:species]), organisms[:,:ncbi_tax_id]))

# load mechanical properties 🕸

mechanical_properties = CSV.read(datadir("silkome/mechanical_properties.csv"), DataFrame)
mechanical_properties = mechanical_properties[:, [ncbi_tax_id, :sex, :country, :toughness, :supercontraction]]

# load proteins 🧬

open(FASTA.Reader, datadir("silkome/spider-silkome-database.v1.prot.fasta")) do reader
    global sequences = [(identifier(record), sequence(record)) for record in reader]
end

getattr(id) = String.(split(id, "|")) |> (t->(family=t[3], genus=t[4], species=t[5], spindroins=t[6:end]))

sequences = DataFrame([(getattr(id)..., seq,) for (id, seq) in sequences])
rename!(sequences, [:family, :genus, :species, :spindroins, :sequence])

sequences[:,:ncbi_tax_id] = [org2id[Tuple(sp)] for sp in eachrow(sequences[:,[1,2,3]])]

sequences_with_properties = innerjoin(sequences, mechanical_properties, on=:ncbi_tax_id, matchmissing=:notequal, makeunique=true)

using Plots, BaHDC


organisms

sequences

mechanical_properties




# visualizing & summarizing the data
# ----------------------------------

# hier kunnen jullie werken aan je case van silkome

histogram(mechanical_properties[:,:toughness])

histogram(mechanical_properties[:,:supercontraction])

describe(mechanical_properties)


# HDC experiments
# ---------------

amino_acid_hdvs = Dict(aa=>hdv() for aa in 'A':'Z')

encode_sequence("AABC", 2, 10_000, amino_acid_hdvs)

# encode all sequences per species

# bin the material properties in 2-3 bins

# use svd to plot in 2 dimensions


X = randn(100, 10_000)  # this is a desciptor matrix

U, s, V = svd(X)

# fraction of the variance explained in 2d
f2d = (s[1] + s[2]) / sum(s)

plot(plot(s, title="spectrum"), scatter(U[:,1], U[:,2]))

