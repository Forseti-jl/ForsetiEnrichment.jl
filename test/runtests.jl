using Test
using ForsetiCore
using ForsetiOmicsCore
using ForsetiEnrichment
using DataFrames
using Statistics
using Random
using Distributions

@testset "ForsetiEnrichment" begin

    @testset "enrich_ora: exact hypergeometric, hand-derived via binomial()" begin
        background = ["g$i" for i in 1:20]
        gene_sets = Dict("setA" => ["g$i" for i in 1:10], "setB" => ["g$i" for i in 11:20])
        sig_genes = vcat(["g$i" for i in 1:6], ["g11", "g12"])  # 6 in setA, 2 in setB

        fit = enrich_ora(sig_genes, gene_sets; background = background)
        @test fit isa OraFit
        @test fit isa ForsetiFit
        @test fit.background_size == 20
        @test fit.n_sig == 8

        t = tidy(fit)
        rowA = t[t.set .== "setA", :]
        @test rowA.k[1] == 6 && rowA.n[1] == 10 && rowA.K[1] == 8 && rowA.N[1] == 20

        # independent combinatorial derivation: P(X >= 6), X ~ Hypergeometric(N=20,K=8,n=10)
        expected_pA = sum(binomial(8, kk) * binomial(12, 10 - kk) for kk in 6:8) / binomial(20, 10)
        @test rowA.p_value[1] ≈ expected_pA atol = 1e-12

        rowB = t[t.set .== "setB", :]
        @test rowB.k[1] == 2
        expected_pB = 1 - (binomial(8, 0) * binomial(12, 10) + binomial(8, 1) * binomial(12, 9)) / binomial(20, 10)
        @test rowB.p_value[1] ≈ expected_pB atol = 1e-12

        # setA is strongly enriched, setB is not
        @test rowA.p_value[1] < 0.1
        @test rowB.p_value[1] > 0.9

        @test issorted(t.p_value)
        g = glance(fit)
        @test g.n_sets[1] == 2
        @test g.background_size[1] == 20
    end

    @testset "enrich_ora: default background, DataFrame/pipe form, OmicsExperiment form" begin
        gene_sets = Dict("setA" => ["g$i" for i in 1:10], "setB" => ["g$i" for i in 11:20])
        sig_genes = vcat(["g$i" for i in 1:6], ["g11", "g12"])

        # no explicit background -> union of all gene_sets values (still 20 genes here)
        fit_default = enrich_ora(sig_genes, gene_sets)
        @test fit_default.background_size == 20

        sig_df = DataFrame(feature = sig_genes)
        fit_pipe = sig_df |> enrich_ora(:feature, gene_sets)
        @test sort(tidy(fit_pipe).p_value) ≈ sort(tidy(fit_default).p_value) atol = 1e-12

        oe = OmicsExperiment(zeros(20, 1), ["g$i" for i in 1:20], ["s1"])
        fit_oe = enrich_ora(oe, sig_genes, gene_sets)
        @test fit_oe.background_size == 20
        @test sort(tidy(fit_oe).p_value) ≈ sort(tidy(fit_default).p_value) atol = 1e-12

        @test_throws ArgumentError enrich_ora(String[], gene_sets)  # no sig genes
    end

    @testset "enrich_gsea: exact enrichment score, hand-derived" begin
        genes = ["g$i" for i in 1:10]
        scores = Float64.(10:-1:1)  # already rank-ordered: g1 highest, g10 lowest

        # hits = top 3 ranks -> perfect early enrichment -> ES = +1.0 exactly
        fit_top = enrich_gsea(genes, scores, Dict("top3" => ["g1", "g2", "g3"]);
                               n_perm = 50, rng = MersenneTwister(1))
        @test fit_top isa GseaFit
        @test fit_top isa ForsetiFit
        @test fit_top.table.es[1] ≈ 1.0 atol = 1e-10
        @test fit_top.table.nes[1] > 0

        # hits = bottom 3 ranks -> perfect late (anti-)enrichment -> ES = -1.0 exactly
        fit_bottom = enrich_gsea(genes, scores, Dict("bottom3" => ["g8", "g9", "g10"]);
                                  n_perm = 50, rng = MersenneTwister(1))
        @test fit_bottom.table.es[1] ≈ -1.0 atol = 1e-10
        @test fit_bottom.table.nes[1] < 0
    end

    @testset "enrich_gsea: permutation p-value is smaller for real signal than for scattered hits" begin
        genes = ["g$i" for i in 1:10]
        scores = Float64.(10:-1:1)
        gene_sets = Dict("top3" => ["g1", "g2", "g3"], "scattered" => ["g2", "g5", "g7", "g9"])

        fit = enrich_gsea(genes, scores, gene_sets; n_perm = 500, rng = MersenneTwister(42))
        t = tidy(fit)
        p_top = t[t.set .== "top3", :p_value][1]
        p_scattered = t[t.set .== "scattered", :p_value][1]
        @test p_top < p_scattered
        @test 0.0 < p_top <= 1.0
        @test 0.0 < p_scattered <= 1.0
        @test issorted(t.p_value)
    end

    @testset "enrich_gsea: reproducible with a seeded rng" begin
        genes = ["g$i" for i in 1:10]
        scores = Float64.(10:-1:1)
        gene_sets = Dict("top3" => ["g1", "g2", "g3"])

        fit1 = enrich_gsea(genes, scores, gene_sets; n_perm = 200, rng = MersenneTwister(7))
        fit2 = enrich_gsea(genes, scores, gene_sets; n_perm = 200, rng = MersenneTwister(7))
        @test fit1.table.p_value == fit2.table.p_value
        @test fit1.table.nes == fit2.table.nes
    end

    @testset "enrich_gsea: DataFrame/pipe form, min_size filtering, errors" begin
        genes = ["g$i" for i in 1:10]
        scores = Float64.(10:-1:1)
        df = DataFrame(feature = genes, score = scores)

        t = df |> enrich_gsea(:feature, :score, Dict("top3" => ["g1", "g2", "g3"]);
                               n_perm = 50, rng = MersenneTwister(1)) |> tidy
        @test t.es[1] ≈ 1.0 atol = 1e-10

        # a gene set with no genes present in the ranked list is skipped entirely...
        @test_throws ArgumentError enrich_gsea(genes, scores, Dict("absent" => ["zzz1", "zzz2"]);
                                                 n_perm = 20, rng = MersenneTwister(1))

        # ...but is fine when at least one other set is usable
        fit_mixed = enrich_gsea(genes, scores,
                                 Dict("absent" => ["zzz1", "zzz2"], "top3" => ["g1", "g2", "g3"]);
                                 n_perm = 20, rng = MersenneTwister(1))
        @test nrow(fit_mixed.table) == 1
        @test fit_mixed.table.set[1] == "top3"
    end

end
