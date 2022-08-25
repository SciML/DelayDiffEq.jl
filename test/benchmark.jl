@testset "I/O" begin
    function check_dataframe(df, wp = nothing)
        # Simple checks (ensures that LaTeX plots are not broken)
        @test names(df) == ["alg", "error", "time"]
        @test nrow(df) == 88
        algs = [
            "RK4",
            "DP5",
            "OwrenZen3",
            "OwrenZen5",
            "Vern6",
            "Vern7",
            "Vern8",
            "Vern9",
            "Rosenbrock23",
            "Rodas4",
            "Rodas5P",
        ]
        @test df.alg == repeat(algs; inner = 8)

        # Check consistency
        if wp !== nothing
            @test wp.names == algs
            for (i, alg) in enumerate(wp.names)
                wp_i = wp[i]
                @test df[(8 * (i - 1) + 1):(8 * i),
                         :] ==
                      DataFrame(; alg = alg,
                                error = wp_i.errors,
                                time = wp_i.times)
            end
        end
    end

    # Check benchmarks in the repo
    @testset "repo" begin
        df = CSV.read(joinpath(@__DIR__, "..",
                               "data",
                               "benchmark.csv"),
                      DataFrame)
        check_dataframe(df)
    end

    # Save benchmark to I/O buffer
    @testset "buffer" begin
        io = IOBuffer()
        wp = benchmark(io)
        df = CSV.read(seekstart(io), DataFrame)
        check_dataframe(df)
    end

    # Save benchmark to file in non-existing path
    @testset "file" begin mktempdir() do dir
        # Desired filename in non-existing path
        file = joinpath(dir, "benchmark",
                        "results.csv")
        @assert !isdir(dirname(file))

        wp = benchmark(file)
        df = CSV.read(file, DataFrame)
        check_dataframe(df)
    end end
end
