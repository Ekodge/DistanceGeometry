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
    @testset "myAtom tests" begin
        include("myAtom_tests.jl")
    end
    @testset "aminoacids tests" begin
        include("aminoacids_tests.jl")
    end
    @testset "Realization tests" begin
        include("realization_tests.jl")
    end
    @testset "Orders tests" begin
        include("orders_tests.jl")
    end
    @testset "DGP tests" begin
        include("dgp_tests.jl")
    end
    @testset "Utils tests" begin
        include("utils_tests.jl")
    end
    @testset "Errorlist tests" begin
        include("errorlist_tests.jl")
    end
end