using Obsea
using Test

@testset "Obsea.jl" begin

    @testset "Structures" begin
        @testset "State" begin
            state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
            @test state.model === 1
            @test state.frequency === 2.0
            @test state.x === 3.0
            @test state.y === 4.0
            @test state.vx === 5.0
            @test state.vy === 6.0
        end
        @testset "Trajectory" begin
            state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
            ∅ = EmptyState()
            trajectory = Trajectory([state, ∅])
            @test trajectory == [state, ∅]
            @test trajectory[1] === state
        end
        @testset "Metadata" begin
            metadata = Metadata(0.1)
            @test metadata.weight === 0.1
        end
        @testset "Particle" begin
            state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
            ∅ = EmptyState()
            trajectory = Trajectory([state, ∅])
            metadata = Metadata(0.1)
            particle = Particle(trajectory, metadata)
            @test particle.trajectory === trajectory
            @test particle.metadata === metadata
            particle.metadata.weight *= 2.0
            @test particle.metadata.weight === 0.2
        end
        @testset "Parameters" begin
            params = Parameters(1.0, 0.1, 0.97, 0.03)
            @test params.T === 1.0
            @test params.q === 0.1
            @test params.ps === 0.97
            @test params.pb === 0.03
        end

        @testset "Grid" begin
            grid = Grid(
                (0.0:1000.0:30000),
                (0.0:5.0:360)[1:end-1],
                range(0.0, 25.0, length=513),
            )
            @test convert(Float64, grid.rrange.step) == 1000.0
            @test grid.arange[end] === 355.0
            @test grid.frange.len === 513
            @test convert(Float64, grid.frange.step) == 25.0 / 512
        end
        @testset "Scan" begin
            scan = Scan(ones(5), ones(4, 3, 2))
            @test scan.yr == ones(5)
            @test scan.ya == ones(4, 3, 2)
        end

    end

    @testset "Motion Model" begin
        @testset "move" begin
            state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
            ∅ = EmptyState()
            params = Parameters(1.0, 0.0, 0.97, 0.03)
            @test move(state, params) == State(1, 2.0, 8.0, 10.0, 5.0, 6.0)
        end
        @testset "logfk" begin
            state = State(1, 2.0, 3.0, 4.0, 5.0, 6.0)
            ∅ = EmptyState()
            params = Parameters(1.0, 0.1, 0.97, 0.03)
            @test logfk(state, state, params) === log(params.ps)  # TODO
            @test logfk(∅, state, params) === log(1.0 - params.ps)
            @test logfk(state, ∅, params) === log(params.pb)
            @test logfk(∅, ∅, params) === log(1.0 - params.pb)
        end
    end

end # testset
