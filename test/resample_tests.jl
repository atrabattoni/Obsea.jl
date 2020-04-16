@testset "resample.jl" begin

    @testset "argsample" begin
        w = ones(10)
        w /= sum(w)
        @test argsample(w, 10) == collect(1:10)
    end

end
