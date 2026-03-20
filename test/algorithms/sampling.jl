import AlgebraicSolving: _real_roots, _sample_points

@testset "Algorithms -> Real roots" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    f = [x + 2 * y + 2 * z - 1, x^2 + 2 * y^2 + 2 * z^2 - x, 2 * x * y + 2 * y * z - y]
    rs = _real_roots(f)
    @test length(rs) == 4
end

@testset "Algorithms -> Univariate sample points" begin
    R, (x,) = polynomial_ring(QQ, ["x"])
    f = [one(R)]
    @test length(_sample_points(f)) == 1
    f = [x + 1, x + 4294967295 // 4294967296, x + 4294967297 // 4294967296]
    @test length(_sample_points(f)) == 4
end

@testset "Algorithms -> Bivariate sample points" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"])
    f = [one(R)]
    @test length(_sample_points(f)) == 1
    f = [x, x + 1]
    @test length(_sample_points(f)) == 3
    f = [y, y + 1]
    @test length(_sample_points(f)) == 3
    f = [x, x + 1, y - x, y + x - 1]
    @test length(_sample_points(f)) >= 9
end

@testset "Algorithms -> Trivariate sample points" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    f = [one(R)]
    @test length(_sample_points(f)) == 1
    f = [x, x + 1]
    @test length(_sample_points(f)) == 3
    f = [y, y + 1]
    @test length(_sample_points(f)) == 3
    f = [x, x + 1, y - x, y + x - 1]
    @test length(_sample_points(f)) >= 9
    f = [x, y, z]
    @test length(_sample_points(f)) == 8
end
