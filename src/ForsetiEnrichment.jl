"""
    ForsetiEnrichment

Gene set enrichment: [`enrich_ora`](@ref) (over-representation /
hypergeometric test) and [`enrich_gsea`](@ref) (a GSEA-style running-sum
enrichment score with permutation-based significance).

```julia
tidy(diffexpr_fit) |> enrich_gsea(:feature, :t, gene_sets) |> tidy
sig_df |> enrich_ora(:feature, gene_sets) |> tidy
```
"""
module ForsetiEnrichment

using DataFrames
using Statistics
using Random
using Distributions
using ForsetiCore
using ForsetiOmicsCore

export OraFit, enrich_ora
export GseaFit, enrich_gsea
export tidy, glance

function _p_adjust_bh(pvalues::AbstractVector{<:Real})
    n = length(pvalues)
    n == 0 && return Float64[]
    order = sortperm(pvalues)
    sorted_p = pvalues[order]
    adjusted = similar(sorted_p, Float64)
    adjusted[n] = clamp(sorted_p[n], 0.0, 1.0)
    for i in (n - 1):-1:1
        adjusted[i] = min(adjusted[i + 1], clamp(sorted_p[i] * n / i, 0.0, 1.0))
    end
    result = similar(adjusted)
    result[order] = adjusted
    return result
end

# ---------------------------------------------------------------------------
# Over-representation analysis (hypergeometric / Fisher's exact test)
# ---------------------------------------------------------------------------

"""
    OraFit <: ForsetiFit

Result of [`enrich_ora`](@ref).
"""
struct OraFit <: ForsetiFit
    table::DataFrame
    background_size::Int
    n_sig::Int
end

"""
    enrich_ora(sig_genes::AbstractVector{<:AbstractString},
               gene_sets::AbstractDict{<:AbstractString,<:AbstractVector{<:AbstractString}};
               background::Union{Nothing,AbstractVector{<:AbstractString}} = nothing) -> OraFit

Over-representation analysis: for each gene set, an exact one-sided
hypergeometric test of whether `sig_genes` overlap it more than expected
by chance, given `background` (the gene universe; if omitted, the union of
all genes across `gene_sets` is used).
"""
function enrich_ora(sig_genes::AbstractVector{<:AbstractString},
                     gene_sets::AbstractDict{<:AbstractString,<:V};
                     background::Union{Nothing,AbstractVector{<:AbstractString}} = nothing) where {V<:AbstractVector{<:AbstractString}}
    bg = background === nothing ? reduce(union, values(gene_sets)) : background
    bg_set = Set(String.(bg))
    sig_set = intersect(Set(String.(sig_genes)), bg_set)
    N = length(bg_set)
    K = length(sig_set)
    K >= 1 || throw(ArgumentError("no significant genes found within the background"))

    rows = NamedTuple[]
    for (name, genes) in gene_sets
        set_in_bg = intersect(Set(String.(genes)), bg_set)
        n = length(set_in_bg)
        n == 0 && continue
        k = length(intersect(set_in_bg, sig_set))
        pval = k == 0 ? 1.0 : ccdf(Hypergeometric(K, N - K, n), k - 1)
        push!(rows, (set = name, k = k, n = n, K = K, N = N, p_value = pval))
    end
    isempty(rows) &&
        throw(ArgumentError("no gene sets had any overlap with the background"))

    table = DataFrame(rows)
    table[!, :adj_p_value] = _p_adjust_bh(table.p_value)
    sort!(table, :p_value)
    return OraFit(table, N, K)
end

"""
    enrich_ora(oe::OmicsExperiment, sig_genes::AbstractVector{<:AbstractString}, gene_sets; kwargs...)

Run [`enrich_ora`](@ref) using every feature of `oe` as the background
gene universe (a common real-world default: only genes actually measured,
not every gene in the genome).
"""
function enrich_ora(oe::OmicsExperiment, sig_genes::AbstractVector{<:AbstractString},
                     gene_sets::AbstractDict; kwargs...)
    return enrich_ora(sig_genes, gene_sets; background = feature_names(oe), kwargs...)
end

"""
    enrich_ora(df::AbstractDataFrame, gene_col::Symbol, gene_sets; kwargs...)

Run [`enrich_ora`](@ref) using `df[!, gene_col]` as the significant gene
list (e.g. a topTable filtered to `adj_p_value < 0.05`).
"""
function enrich_ora(df::AbstractDataFrame, gene_col::Symbol, gene_sets::AbstractDict; kwargs...)
    return enrich_ora(String.(df[!, gene_col]), gene_sets; kwargs...)
end

"""
    enrich_ora(gene_col::Symbol, gene_sets; kwargs...)

Pipe-curried form: `sig_df |> enrich_ora(:feature, gene_sets)`.
"""
enrich_ora(gene_col::Symbol, gene_sets::AbstractDict; kwargs...) =
    df -> enrich_ora(df, gene_col, gene_sets; kwargs...)

ForsetiCore.tidy(fit::OraFit) = fit.table

function ForsetiCore.glance(fit::OraFit; alpha::Real = 0.05)
    return DataFrame(n_sets = [nrow(fit.table)], background_size = [fit.background_size],
                      n_sig_genes = [fit.n_sig], n_enriched = [count(<(alpha), fit.table.adj_p_value)])
end

# ---------------------------------------------------------------------------
# GSEA-style enrichment score with permutation-based significance
# ---------------------------------------------------------------------------

"""
    GseaFit <: ForsetiFit

Result of [`enrich_gsea`](@ref).
"""
struct GseaFit <: ForsetiFit
    table::DataFrame
    n_perm::Int
end

function _enrichment_score(scores::Vector{Float64}, hit::BitVector, weight::Real)
    N = length(scores)
    Nh = count(hit)
    Nh == 0 && return 0.0
    weights = abs.(scores) .^ weight
    hit_sum = sum(weights[hit])
    miss_penalty = 1.0 / (N - Nh)
    running = 0.0
    max_dev = 0.0
    for i in 1:N
        if hit[i]
            running += hit_sum == 0 ? 0.0 : weights[i] / hit_sum
        else
            running -= miss_penalty
        end
        abs(running) > abs(max_dev) && (max_dev = running)
    end
    return max_dev
end

"""
    enrich_gsea(genes::AbstractVector{<:AbstractString}, scores::AbstractVector{<:Real},
                gene_sets::AbstractDict{<:AbstractString,<:AbstractVector{<:AbstractString}};
                n_perm = 1000, min_size = 1, weight = 1.0, rng = Random.default_rng()) -> GseaFit

GSEA-style enrichment: `genes` are ranked by `scores` (descending), and for
each gene set a running-sum enrichment score (ES) is computed (the
classic weighted Kolmogorov-Smirnov statistic from Subramanian et al.
2005, with exponent `weight`), then a p-value is obtained by comparing ES
to a null distribution built from `n_perm` random gene sets of the same
size (gene-set permutation).
"""
function enrich_gsea(genes::AbstractVector{<:AbstractString}, scores::AbstractVector{<:Real},
                      gene_sets::AbstractDict{<:AbstractString,<:V};
                      n_perm::Int = 1000, min_size::Int = 1, weight::Real = 1.0,
                      rng::AbstractRNG = Random.default_rng()) where {V<:AbstractVector{<:AbstractString}}
    ForsetiCore.require_same_length(genes, scores; xname = "genes", yname = "scores")
    order = sortperm(scores; rev = true)
    sorted_genes = String.(genes[order])
    sorted_scores = Float64.(collect(scores[order]))
    N = length(sorted_genes)
    gene_pos = Dict(g => i for (i, g) in enumerate(sorted_genes))

    rows = NamedTuple[]
    for (name, set_genes) in gene_sets
        hit = falses(N)
        for g in set_genes
            haskey(gene_pos, g) && (hit[gene_pos[g]] = true)
        end
        Nh = count(hit)
        Nh < min_size && continue

        es = _enrichment_score(sorted_scores, hit, weight)

        perm_es = Vector{Float64}(undef, n_perm)
        for b in 1:n_perm
            perm_hit = falses(N)
            perm_hit[randperm(rng, N)[1:Nh]] .= true
            perm_es[b] = _enrichment_score(sorted_scores, perm_hit, weight)
        end

        same_sign = es >= 0 ? filter(>=(0), perm_es) : filter(<(0), perm_es)
        nes = isempty(same_sign) ? es : es / mean(abs.(same_sign))
        n_extreme = es >= 0 ? count(>=(es), perm_es) : count(<=(es), perm_es)
        pval = (n_extreme + 1) / (n_perm + 1)

        push!(rows, (set = name, n_genes = Nh, es = es, nes = nes, p_value = pval))
    end
    isempty(rows) &&
        throw(ArgumentError("no gene sets met min_size = $min_size after matching to the ranked gene list"))

    table = DataFrame(rows)
    table[!, :adj_p_value] = _p_adjust_bh(table.p_value)
    sort!(table, :p_value)
    return GseaFit(table, n_perm)
end

"""
    enrich_gsea(df::AbstractDataFrame, gene_col::Symbol, score_col::Symbol, gene_sets; kwargs...)

Run [`enrich_gsea`](@ref) using `df[!, gene_col]`/`df[!, score_col]` as
the ranked gene list, e.g. straight from a differential-expression
topTable: `tidy(diffexpr_fit) |> enrich_gsea(:feature, :t, gene_sets)`.
"""
function enrich_gsea(df::AbstractDataFrame, gene_col::Symbol, score_col::Symbol,
                      gene_sets::AbstractDict; kwargs...)
    return enrich_gsea(String.(df[!, gene_col]), df[!, score_col], gene_sets; kwargs...)
end

"""
    enrich_gsea(gene_col::Symbol, score_col::Symbol, gene_sets; kwargs...)

Pipe-curried form: `topTable |> enrich_gsea(:feature, :t, gene_sets)`.
"""
enrich_gsea(gene_col::Symbol, score_col::Symbol, gene_sets::AbstractDict; kwargs...) =
    df -> enrich_gsea(df, gene_col, score_col, gene_sets; kwargs...)

ForsetiCore.tidy(fit::GseaFit) = fit.table

function ForsetiCore.glance(fit::GseaFit; alpha::Real = 0.05)
    return DataFrame(n_sets = [nrow(fit.table)], n_perm = [fit.n_perm],
                      n_sig = [count(<(alpha), fit.table.adj_p_value)])
end

end # module ForsetiEnrichment
