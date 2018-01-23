using DelayDiffEq, DiffEqProblemLibrary, Base.Test

## Special test

# check constant extrapolation with problem with vanishing delays at t = 0
prob = DDEProblem((du,u,h,p,t) -> -h(t/2)[1], t -> [1.0], [1.0], (0.0, 10.0), nothing, [])
solve(prob, MethodOfSteps(RK4()))

## General tests

# history function
function h(t, idxs=nothing)
    h(t,Val{0},idxs)
end
function h(val::AbstractArray, t, idxs=nothing)
    h(val,t,Val{0},idxs)
end
function h(t, ::Type{Val{0}}, idxs=nothing)
    if typeof(idxs) <: Void
        return [t; -t]
    else
        return [t; -t][idxs]
    end
end
function h(val::AbstractArray, t, ::Type{Val{0}}, idxs=nothing)
    if typeof(idxs) <: Void
        val[1] = t
        val[2] = -t
    else
        val .= [t; -t][idxs]
    end
end

# ODE integrator
prob = ODEProblem(DiffEqProblemLibrary.f_2dlinear, ones(2), (0.0, 1.0))
integrator = init(prob, Tsit5())

# combined history function
history = DelayDiffEq.HistoryFunction(h, integrator.sol, integrator)

@which h(-0.5, nothing)

# test evaluation of history function
for idxs in (nothing, [2])
    @test history(-0.5, Val{0}, idxs) == h(-0.5, idxs)

    val = idxs == nothing ? zeros(2) : [0.0]
    val2 = deepcopy(val)

    history(val, -0.5, Val{0}, idxs)
    h(val2, -0.5, idxs)

    @test val == val2
end

# test constant extrapolation
for (deriv, idxs) in Iterators.product((Val{0}, Val{1}), (nothing, [2]))
    integrator.isout = false

    extrapolation = deriv == Val{0} ?
        (idxs == nothing ? integrator.u : integrator.u[[2]]) :
        (idxs == nothing ? zeros(2) : [0.0])

    @test history(0.5, deriv, idxs) == extrapolation && integrator.isout

    integrator.isout = false
    @test history(nothing, 0.5, deriv, idxs) == extrapolation && integrator.isout

    integrator.isout = false
    val = 1 .- extrapolation # ensures that val ≠ extrapolation
    history(val, 0.5, deriv, idxs)

    @test val == extrapolation && integrator.isout
end

# add step to integrator
OrdinaryDiffEq.loopheader!(integrator)
OrdinaryDiffEq.perform_step!(integrator, integrator.cache)
integrator.t = integrator.dt
@test 0.01 < integrator.t < 1
@test integrator.sol.t[end] == 0

# test integrator interpolation
for (deriv, idxs) in Iterators.product((Val{0}, Val{1}), (nothing, [2]))
    integrator.isout = false

    @test history(0.01, deriv, idxs) ==
        OrdinaryDiffEq.current_interpolant(0.01, integrator, idxs, deriv) &&
        integrator.isout

    integrator.isout = false
    val = idxs == nothing ? zeros(2) : [0.0]
    val2 = deepcopy(val)
    history(val, 0.01, deriv, idxs)
    OrdinaryDiffEq.current_interpolant!(val2, 0.01, integrator, idxs, deriv)

    @test val == val2 && integrator.isout
end

# add step to solution
integrator.t = 0
OrdinaryDiffEq.loopfooter!(integrator)
@test integrator.t == integrator.sol.t[end]

# test solution interpolation
for (deriv, idxs) in Iterators.product((Val{0}, Val{1}), (nothing, [2]))
    @test history(0.01, deriv, idxs) ==
        integrator.sol.interp(0.01, idxs, deriv, integrator.p) &&
        !integrator.isout

    val = idxs == nothing ? zeros(2) : [0.0]
    val2 = deepcopy(val)
    history(val, 0.01, deriv, idxs)
    integrator.sol.interp(val2, 0.01, idxs, deriv, integrator.p)

    @test val == val2 && !integrator.isout
end

# test integrator extrapolation
for (deriv, idxs) in Iterators.product((Val{0}, Val{1}), (nothing, [2]))
    integrator.isout = false

    @test history(1, deriv, idxs) ==
        OrdinaryDiffEq.current_interpolant(1, integrator, idxs, deriv) &&
        integrator.isout

    integrator.isout = false
    val = idxs == nothing ? zeros(2) : [0.0]
    val2 = deepcopy(val)
    history(val, 1, deriv, idxs)
    OrdinaryDiffEq.current_interpolant!(val2, 1, integrator, idxs, deriv)

    @test val == val2 && integrator.isout
end
