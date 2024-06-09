using Test
include("utilstest.jl")

@testset "DGP tests" begin
    @testset "distance tests" begin
        include("distance_tests.jl")
    end
    @testset "dimensionality tests" begin
        include("dimensionality_tests.jl")
    end
    @testset "distlist tests" begin
        include("distlist_tests.jl")
    end
    @testset "aminoacids tests" begin
        include("aminoacids_tests.jl")
    end
    @testset "myAtom tests" begin
        include("myAtom_tests.jl")
    end
end