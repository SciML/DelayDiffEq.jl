@testset "solution" begin
    # Compute the solution
    sol = mackey_glass()

    # Check success and desired properties
    @test sol.retcode === :Success
    @test extrema(sol.t) == (0.0, 600.0)
end

@testset "I/O" begin
    function check_dataframe(df, sol)
        # Compute the desired solution
        testsol = mackey_glass()

        # Check that the data frame is correct
        @test names(df) == ["t", "x", "xtau"]
        @test df.t == 300:0.1:600
        @test df.x == Vector(testsol(300:0.1:600))
        @test df.xtau ==
              Vector(testsol(298:0.1:598))

        # Check consistency
        if sol !== nothing
            @test sol.retcode === :Success
            @test sol.t == testsol.t
            @test sol.u == testsol.u
            @test typeof(sol.prob) ===
                  typeof(testsol.prob)
            @test DataFrame(sol) ==
                  DataFrame(testsol)
        end
    end

    # Save solution to I/O buffer
    @testset "buffer" begin
        io = IOBuffer()
        sol = mackey_glass(io)
        df = CSV.read(seekstart(io), DataFrame)
        check_dataframe(df, sol)
    end

    # Save solution to file in non-existing path
    @testset "file" begin mktempdir() do dir
        # Desired filename in non-existing path
        file = joinpath(dir, "mackey_glass",
                        "sol.csv")
        @assert !isdir(dirname(file))

        sol = mackey_glass(file)
        df = CSV.read(file, DataFrame)
        check_dataframe(df, sol)
    end end
end
