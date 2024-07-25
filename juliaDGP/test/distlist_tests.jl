struct TestDistList{T} <: DistList{T}
    K::Int64
    distlist::Dict{Tuple{T,T},Distance}
end

# Define some test data
distances_dict = Dict{Tuple{Int,Int}, Distance}(
    (1, 2) => Distance(0.0, 1.0),
    (1, 3) => Distance(0.0, 2.0),
    (2, 3) => Distance(0.0, 3.0)
)

distances_dict2 = Dict{Tuple{Int,Int}, Distance}(
    (1, 2) => Distance(1.0, 1.0),
    (1, 3) => Distance(2.0, 2.0),
    (2, 3) => Distance(0.0, 3.0),
    (3, 4) => Distance(1.0,5.0)
)

@testset "distlist functions" begin
    @testset "nd" begin
        dl = TestDistList(3, distances_dict)
        @test nd(dl) == 3
    end

    @testset "nexact" begin
        dl = TestDistList(3, distances_dict)
        @test nexact(dl) == 0
        dl = TestDistList(3, distances_dict2)
        @test nexact(dl) == 2
    end

    @testset "ninterval" begin
        dl = TestDistList(3, distances_dict)
        @test ninterval(dl) == 3
        dl = TestDistList(3, distances_dict2)
        @test ninterval(dl) == 2
    end

    @testset "elements" begin
        dl = TestDistList(3, distances_dict)
        @test elements(dl) == [1, 2, 3]
        dl = TestDistList(3, distances_dict2)
        @test elements(dl) == [1, 2, 3, 4]
    end

    @testset "nelements" begin
        dl = TestDistList(3, distances_dict)
        @test nelements(dl) == 3
        dl = TestDistList(3, distances_dict2)
        @test nelements(dl) == 4
    end

    @testset "details" begin
    dl = TestDistList(3, distances_dict)
    @test print_to_string_stdout(details,dl) == "(1,2) => 0.5\n(1,3) => 1.0\n(2,3) => 1.5\n"
    end

    @testset "sublist" begin
        dl = TestDistList(3, distances_dict)
        @test sublist(dl, [1, 2]) == [Distance(0.0, 1.0)]
        @test sublist(dl, [1, 2, 3]) == [Distance(0.0, 1.0), Distance(0.0, 2.0), Distance(0.0, 3.0)]
        # Test error cases :
        # Error : Elements of the given subset are not in the Distlist instance
        @test_throws ArgumentError sublist(dl, [1, 2, 3, 4])
    end

    @testset "is_clique" begin
        dl = TestDistList(3, distances_dict)
        @test is_clique(dl, [1, 2, 3]) == true
        @test is_clique(dl, [1, 2]) == true
        @test is_clique(dl, [1]) == true

        dl = TestDistList(3, distances_dict2)
        @test is_clique(dl, [1, 2, 3, 4]) == false
        @test is_clique(dl, [1, 2, 4]) == false
        @test is_clique(dl, [1, 2]) == true
        @test is_clique(dl, [1]) == true
    end

    @testset "graph" begin
        dl = TestDistList(3, distances_dict)
        g = graph(dl)
        @test NE(g) == 3
        @test NV(g) == 3
        @test vlist(g) == [1, 2, 3]
        @test elist(g) == [(1, 2), (1, 3), (2, 3)]
    end
end