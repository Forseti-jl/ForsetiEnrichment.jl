# ForsetiEnrichment.jl

Gene set enrichment and over-representation analysis.

Part of the [Forseti](https://github.com/Forseti-jl) statistical analysis package family.

## Usage

```julia
using ForsetiCore, ForsetiOmicsCore, ForsetiEnrichment, DataFrames

gene_sets = Dict("pathwayA" => ["g1", "g2", "g3"], "pathwayB" => ["g4", "g5"])

# over-representation analysis: exact hypergeometric test
sig_df = DataFrame(feature = ["g1", "g2", "g3"])
sig_df |> enrich_ora(:feature, gene_sets) |> tidy      # set, k, n, K, N, p_value, adj_p_value
enrich_ora(oe, sig_genes, gene_sets)                    # background = every feature of oe

# GSEA-style: rank genes by a score (e.g. straight from a diff_expr topTable),
# running-sum enrichment score + permutation-based p-value
tidy(diffexpr_fit) |> enrich_gsea(:feature, :t, gene_sets) |> tidy  # set, n_genes, es, nes, p_value, adj_p_value
```

`enrich_gsea` computes the classic weighted running-sum enrichment score
(Subramanian et al. 2005) exactly, and gets significance from `n_perm`
random gene-set permutations rather than the original paper's asymptotic
corrections — pass a seeded `rng` for reproducibility.

## API Reference

```@autodocs
Modules = [ForsetiEnrichment]
```
