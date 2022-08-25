@testset "prediction" begin
    p = [2.2, 1.0, 2.0, 0.4]
    x = DelayDiffEqPaper.predict_lv(p)
    @test size(x) == (2, 51)
end

@testset "Jacobian" begin
    # Compute Jacobian of Lotka-Volterra model
    (; forwarddiff, zygote) = jacobian_lv()

    # Check dimensions of Jacobians
    @test size(forwarddiff) == (102, 4)
    @test size(zygote) == (102, 4)

    # Check that forward-mode and reverse-mode match
    @test forwarddiff ≈ zygote

    # Compute Jacobian with finite-differencing method
    p = [2.2, 1.0, 2.0, 0.4]
    J = FiniteDiff.finite_difference_jacobian(DelayDiffEqPaper.predict_lv,
                                              p)

    # Compare with ForwardDiff and Zygote
    @test forwarddiff≈J rtol=1e-4
    @test zygote≈J rtol=1e-4

    # At the initial time point, derivatives of both
    # states with respect to all parameters should be zero
    @test all(iszero, forwarddiff[1:2, :])
    @test all(iszero, zygote[1:2, :])
end
