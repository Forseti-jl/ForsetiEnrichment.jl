# ForsetiEnrichment.jl

Gene set enrichment and over-representation analysis.

Part of the [Forseti](https://github.com/natapol) statistical analysis package family.

## Dependencies

This package depends on the following sibling Forseti packages, which are not
yet registered and must be added via local dev paths:

- `ForsetiCore`
- `ForsetiOmicsCore`

## Local development

```julia
using Pkg
Pkg.develop([PackageSpec(path="../ForsetiCore.jl"), PackageSpec(path="../ForsetiOmicsCore.jl")])
Pkg.instantiate()
Pkg.test()
```
