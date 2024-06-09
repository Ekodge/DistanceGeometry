struct dimensionalityTest <: Dimensionality
    K::Int64;
    dimensionalityTest(k::Int64) = new(k)
    end

struct dimensionalityTest2 <: Dimensionality
    K::Int64;
    dimensionalityTest2(k::Int64) = new(k)
    end

@testset "dimension" begin
    d = dimensionalityTest(4)
    @test dimension(d) == 4
    d = dimensionalityTest2(5)
    @test dimension(d) == 5
    end

@testset "are_compatible" begin
    d1 = dimensionalityTest(4)
    d2 = dimensionalityTest(4)
    @test are_compatible(d1, d2) == true
    d2 = dimensionalityTest(5)
    @test are_compatible(d1, d2) == false
    d3 = dimensionalityTest2(5)
    @test are_compatible(d1, d3) == false
    @test are_compatible(d2, d3) == true
    end