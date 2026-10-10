using Test
using LinearAlgebra
using Statistics
using Random
using DataFrames
import BioStructures
using RibosomeAnalysis

include("structures/test_structures.jl")
include("binding/test_binding.jl")
include("interactions/test_interactions.jl")
include("geometry/test_geometry.jl")
include("alignment/test_alignment.jl")
include("alignment/test_comparison.jl")
include("dynamics/test_dynamics.jl")
include("statistics/test_statistics.jl")
include("io/test_results.jl")
include("validation/test_mdanalysis_reference.jl")
include("quality/test_quality.jl")
