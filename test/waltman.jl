@testset "solution" begin
    # Compute the solution
    sol = waltman()

    # Check success and desired properties
    @test sol.retcode === :Success
    @test sol.t == 0:0.1:300
    @test typeof(sol.prob) ===
          typeof(DiffEqBase.get_concrete_problem(DDEProblemLibrary.prob_dde_RADAR5_waltman_5,
                                                 true))

    # Load the solution computed with RADAR5
    df = CSV.read(joinpath(@__DIR__, "radar5",
                           "waltman.csv"),
                  DataFrame)
    @assert df.timestamp == sol.t
    @assert names(df) == vcat("timestamp",
                 ["value$i" for i in 1:6])

    # Compare different errors
    sol_u = reduce(hcat, sol.u)
    sol_u_14 = sol_u[1:4, :]
    sol_u_56 = sol_u[5:6, :]
    df_u_14 = Matrix(df[:,
                        ["value$i" for i in 1:4]])'
    df_u_56 = Matrix(df[:, ["value5", "value6"]])'

    # Final errors
    @test mean(abs,
               sol_u_14[:, end] - df_u_14[:, end]) <
          3e-10
    @test mean(abs,
               sol_u_56[:, end] - df_u_56[:, end]) <
          3e-4

    # l∞ errors
    @test maximum(abs.(sol_u_14 .- df_u_14)) < 2e-9
    @test maximum(abs.(sol_u_56 .- df_u_56)) < 2e-1

    # l2 errors
    @test sqrt(mean(abs2.(sol_u_14 .- df_u_14))) <
          3e-10
    @test sqrt(mean(abs2.(sol_u_56 .- df_u_56))) <
          3e-3
end

@testset "I/O" begin
    function check_dataframe(df, sol)
        # Compute the desired solution
        testsol = waltman()

        # Check that the data frame is correct
        @test df == DataFrame(testsol)

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
        sol = waltman(io)
        df = CSV.read(seekstart(io), DataFrame)
        check_dataframe(df, sol)
    end

    # Save solution to file in non-existing path
    @testset "file" begin mktempdir() do dir
        # Desired filename in non-existing path
        file = joinpath(dir, "waltman", "sol.csv")
        @assert !isdir(dirname(file))

        sol = waltman(file)
        df = CSV.read(file, DataFrame)
        check_dataframe(df, sol)
    end end
end
