@testset "Distance constructors" begin

    @testset "Generic constructor" begin
        d = Distance(1.0, 2.0, 1)
        @test d.lb == 1.0
        @test d.ub == 2.0
        @test d.uncertainty == 1
        # Test error cases :
        # Error : Lower bound is larger than upper bound
        @test_throws ArgumentError Distance(2.0, 1.0, 1)
        # Error : The uncertainty level cannot be negative
        @test_throws ArgumentError Distance(1.0, 2.0, -1)
    end

    @testset "Generic constructor with boolean" begin
        d = Distance(1.0, 2.0, true)
        @test d.lb == 1.0
        @test d.ub == 2.0
        @test d.uncertainty == 1
        d = Distance(4.0, 6.0, false)
        @test d.lb == 4.0
        @test d.ub == 6.0
        @test d.uncertainty == 0
        # Test error cases :
        # Error : Lower bound is larger than upper bound
        @test_throws ArgumentError Distance(2.0, 1.0,true)
    end

    @testset "Constructor without uncertainty" begin
        d = Distance(1.0, 2.0)
        @test d.lb == 1.0
        @test d.ub == 2.0
        @test d.uncertainty == 0
        # Test error cases :
        # Error : Lower bound is larger than upper bound
        @test_throws ArgumentError Distance(2.0, 1.0)
    end

    @testset "Constructor for exact distances" begin
        d = Distance(1.0)
        @test d.lb == 1.0
        @test d.ub == 1.0
        @test d.uncertainty == 0
    end

    @testset "Constructor from a pair of Vector types" begin
        u = [1.0, 2.0]
        v = [4.0, 6.0]
        d = Distance(u, v)
        #sqrt((1-4)^2 + (2-6)^2) = sqrt(9 + 16) = sqrt(25) = 5)
        @test d.lb ≈ 5
        @test d.ub ≈ 5
        # Test error cases :
        # Error : The two Vector instances need to have the same length
        @test_throws ArgumentError Distance([1.0, 2.0], [3.0])
    end
end

@testset "Distance type various overrides" begin

    @testset "Override of the equality operator" begin
        d1 = Distance(1.0, 2.0, 1)
        d2 = Distance(1.0, 2.0, 1)
        @test d1 == d2
        d2 = Distance(1.0, 3.0, 1)
        # d1.lb == d2.lb but d1.ub != d2.ub
        @test d1 != d2
        d2 = Distance(2.0, 2.0, 1)
        # d1.ub == d2.ub but d1.lb != d2.lb
        @test d1 != d2
        d2 = Distance(1.0, 2.0, 0)
        # d1.uncertainty != d2.uncertainty
        @test d1 != d2
    end

    @testset "Override of the hash function" begin
        d1 = Distance(1.0, 2.0, 0)
        @test hash(d1) == hash(1.0) + hash(2.0)
        @test hash(d1) != (hash(1.0) + hash(2.0)) << 1
        d1 = Distance(1.0, 2.0, 1)
        @test hash(d1) == (hash(1.0) + hash(2.0)) << 1
        @test hash(d1) != hash(1.0) + hash(2.0)
        d2 = Distance(1.0, 2.0, 1)
        @test hash(d1) == hash(d2)
    end

    @testset "Override of the show function" begin
        print_to_string(f::Function, d::Distance) = (io = IOBuffer(); f(io,d); String(take!(io)))
        # Test the show function
        d = Distance(1.0, 1.0, 0)
        @test print_to_string(show,d) == "Distance(1.0)"
        d = Distance(1.0, 2.0, 0)
        @test print_to_string(show,d)  == "Distance(1.0,2.0)"
        d = Distance(1.0, 2.0, 3)
        @test print_to_string(show,d)  == "Distance(1.0,2.0,uncertainty level 3)"
    end
end

@testset "Distance Functions" begin

    @testset "is_exact" begin
        d = Distance(1.0)
        @test is_exact(d) == true
        d = Distance(1.0, 2.0, 0)
        @test is_exact(d) == false
    end

    @testset "are_exact" begin
        d1 = Distance(1.0)
        d2 = Distance(1.0)
        @test are_exact([d1, d2]) == true
        d2 = Distance(1.0, 2.0, 0)
        @test are_exact([d1, d2]) == false
    end

    @testset "nexact" begin
        d1 = Distance(1.0)
        d2 = Distance(1.0, 1.0, 0)
        d3 = Distance(1.0, 2.0, 1)
        @test nexact([d1, d2, d3]) == 2
        @test nexact([d1, d2]) == 2
        @test nexact([d1, d3]) == 1
    end

    @testset "is_interval" begin
        d = Distance(1.0, 2.0, 1)
        @test is_interval(d) == true
        d = Distance(1.0)
        @test is_interval(d) == false
    end

    @testset "are_intervals" begin
        d1 = Distance(1.0, 2.0, 1)
        d2 = Distance(1.0, 2.0, 0)
        @test are_intervals([d1]) == true
        @test are_intervals([d2]) == true
        @test are_intervals([d1, d2]) == true
        d2 = Distance(1.0)
        @test are_intervals([d2]) == false
        @test are_intervals([d1, d2]) == false
    end

    @testset "ninterval" begin
        d1 = Distance(1.0, 2.0, 1)
        d2 = Distance(1.0, 2.0, 0)
        d3 = Distance(1.0)
        @test ninterval([d1, d2, d3]) == 2
        @test ninterval([d1, d2]) == 2
        @test ninterval([d1]) == 1
    end

    @testset "is_uncertain" begin
        d = Distance(1.0, 2.0, 1)
        @test is_uncertain(d) == true
        d = Distance(1.0)
        @test is_uncertain(d) == false
    end

    @testset "are_uncertain" begin
        d1 = Distance(1.0, 2.0, 1)
        d2 = Distance(1.0, 2.0, 0)
        @test are_uncertain([d1]) == true
        @test are_uncertain([d2]) == false
        @test are_uncertain([d1, d2]) == false
        d2 = Distance(4.0, 6.0, 1)
        @test are_uncertain([d2]) == true
        @test are_uncertain([d1, d2]) == true
    end

    @testset "nuncertain" begin
        d1 = Distance(1.0, 2.0, 1)
        d2 = Distance(1.0, 2.0, 0)
        d3 = Distance(4.0, 6.0, 1)
        @test nuncertain([d1, d2, d3]) == 2
        @test nuncertain([d1, d3]) == 2
        @test nuncertain([d1]) == 1
        @test nuncertain([d1, d2]) == 1
    end

    @testset "range" begin
        d = Distance(1.0, 2.0, 1)
        @test range(d) == 1.0
        d = Distance(4.0, 6.0, 0)
        @test range(d) == 2.0
        d = Distance(1.0, 2.0, true)
        @test range(d) == 1.0
        d = Distance(1.0)
        @test range(d) == 0.0
    end

    @testset "value" begin
        d = Distance(1.0, 2.0, 1)
        @test value(d) == 1.5
        d = Distance(4.0, 6.0, 0)
        @test value(d) == 5.0
        d = Distance(1.0, 2.0, true)
        @test value(d) == 1.5
        d = Distance(1.0)
        @test value(d) == 1.0
    end

    @testset "intersect" begin
        d1 = Distance(1.0)
        d2 = Distance(1.0)
        # d1 and d2 are exact distances with same lb
        @test intersect(d1,d2) == Distance(1.0, 1.0)
        # d1 is exact and d2 is an interval and lb(d2) <= lb(d1) && lb(d1) <= ub(d2)
        d2 = Distance(1.0, 2.0)
        @test intersect(d1,d2) == Distance(1.0, 1.0)
        # the other way around
        @test intersect(d2,d1) == Distance(1.0, 1.0)
        # d1 and d2 are intervals and lb(d2) <= lb(d1) && lb(d1) <= ub(d2)
        d1 = Distance(1.5, 3.0)
        @test intersect(d1,d2) == Distance(1.5, 2.0)
        # same but no intersection
        d1 = Distance(4.0, 6.0)
        @test intersect(d1,d2) == nothing
        # Test error cases :
        # Error : It is not possible to intersect uncertain distances
        @test_throws ArgumentError intersect(Distance(1.0, 2.0, 1), Distance(1.0, 3.0, 1))
        @test_throws ArgumentError intersect(Distance(1.0, 2.0), Distance(1.0, 3.0, 1))
    end

    @testset "distance" begin
        d1 = Distance(1.0, 2.0)
        d2 = Distance(1.5, 3.0)
        # intersect != nothing
        @test distance(d1, d2) == 0.0
        # return min(abs(1.0 - 6.0),abs(4.0 - 2.0))
        d2 = Distance(4.0, 6.0)
        @test distance(d1, d2) != 5.0
        @test distance(d1, d2) == 2.0
    end

    @testset "with_higher_uncertainty" begin
        d = Distance(1.0, 2.0, 1)
        @test with_higher_uncertainty(d) == Distance(1.0, 2.0, 2)
        # the hash code is the same of the original instance
        #@test hash((with_higher_uncertainty(d))) = hash(d)
        # Test error cases :
        # Error : This Distance type does not hold an uncertain distance
        @test_throws ArgumentError with_higher_uncertainty(Distance(1.0, 2.0))
    end
end