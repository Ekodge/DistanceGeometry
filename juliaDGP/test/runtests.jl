using Test

@testset "DGP tests" begin
    @testset "Distance tests" begin
        include("distance_tests.jl")
    end
end