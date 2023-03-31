module BaHDC

# a module for working with HDC
# hier is jullie source


export hdv, bundle, shift, cossim, hamming, encode_sequence


hdv(;n=10_000) = rand((-1, 1), n)
hdv(m; n=10_000) = rand((-1, 1), m, n)

bundle(x::AbstractVector{Int}, y::AbstractVector{Int}) = sign.(x .+ y)

bundle(xs::AbstractVector{Int}...) = sign.(reduce(.+, xs))

bundle(X::AbstractMatrix{Int}...) = sign.(vec(reduce(+, X, dims=1)))

shift(x, k=1) = circshift(x, k)

cossim(x, y) = dot(x, y) / (norm(x) * norm(y))

# Hamming afstand

hamming(x, y) = sum(x .!= y)


function encode_sequence(sequence, k, N, hdv_dict)
    n = length(sequence)
    encoding = zeros(Int, N)
    temp = similar(encoding)
    for i in 1:n-k
        temp .= hdv_dict[sequence[i]]
        for j in 1:k-1
            v = hdv_dict[sequence[i+j]]
            for pos in 1:N-j
                @inbounds temp[pos] *= v[pos+j]
            end
        end
        encoding .+= temp
    end
    encoding .= sign.(encoding)
    return encoding
end

using LinearAlgebra


        

end