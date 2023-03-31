#=
Created on 30/03/2023 12:25:56
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Loader for the rough silk data loader.

Data from 

K. Arakawa et al., “1000 spider silkomes: Linking sequences to silk physical properties,” Science Advances, vol. 8, no. 41, p. eabo6043, Oct. 2022, doi: 10.1126/sciadv.abo6043.

=#

using DrWatson
@quickactivate "BaHDC"
using DataFrames, FASTX, CSV

# load organsims 🕷

organisms = CSV.read(datadir("silkome/organisms.csv"), DataFrame)

# load mechanical properties 🕸

mechanical_properties = CSV.read(datadir("silkome/mechanical_properties.csv"), DataFrame)
mechanical_properties = mechanical_properties[:, [:family, :genus, :species, :sex, :country, :toughness, :supercontraction]]

# load proteins 🧬

open(FASTA.Reader, datadir("silkome/spider-silkome-database.v1.prot.fasta")) do reader
    global sequences = [(identifier(record), sequence(record)) for record in reader]
end

getattr(id) = String.(split(id, "|")) |> (t->(family=t[3], genus=t[4], species=t[5], spindroins=t[6:end]))

sequences = DataFrame([(getattr(id)..., seq,) for (id, seq) in sequences] )

