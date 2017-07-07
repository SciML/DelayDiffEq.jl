using DelayDiffEq, DiffEqBase, OrdinaryDiffEq, Base.Test

lags = [1]
f = function (t,u,h)
    - h(t-1)
end
h = (t) -> 0.0

function (p::typeof(f))(::Type{Val{:analytic}}, t, u0)
    t = t-1
    if t<0
        return u0
    elseif t<1
        return u0 - t*u0
    elseif t<2
        return 0.5*(3u0 -4*t*u0 + t^2 * u0)
    elseif t<3
        return 1/6*(17u0 -24t*u0 + 9t^2*u0 -t^3 * u0)
    elseif t<4
        return 1/24*(149u0 - 204*t*u0 + 90t^2*u0 - 16t^3*u0 + t^4*u0)
    elseif t<5
        return 1/120*(1769u0 - 2300t*u0 + 1090t^2*u0 - 240t^3*u0 + 25t^4*u0 - t^5*u0)
    elseif t<6
        return 1/720*(26239u0 - 32550t*u0 + 15915t^2*u0 - 3940t^3*u0 + 525t^4*u0 - 36t^5*u0 + t^6*u0)
    elseif t<7
        return (463609u0 - 554442t*u0 + 274701t^2*u0 - 72940t^3*u0 + 11235t^4*u0 - 1008t^5*u0 + 49t^6*u0 - t^7*u0)/5040
    elseif t<8
        return (9473673u0 - 11023880t*u0 + 5491780t^2*u0 - 1524712t^3*u0 + 257950t^4*u0 - 27272t^5*u0 + 1764t^6*u0 - 64t^7*u0 + t^8*u0)/40320
    elseif t<9
        return (219480785u0 - 250209864t*u0 + 124923492t^2*u0 - 35742504t^3*u0 +6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
    elseif t<=10
        return (219480785*u0 - 250209864*t*u0 + 124923492*t^2*u0 - 35742504t^3*u0 + 6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
    else
        error("This analytical solution is only valid on [-infty,11]")
    end
end


prob = ConstantLagDDEProblem(f, h, 1.0, lags, (0.0, 10.0); iip=DiffEqBase.isinplace(f, 4))
alg = MethodOfSteps(BS3(); constrained=true)

dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test maximum(sol.errors[:final]) < 1e-4

h = (t) -> [0.0]
prob = ConstantLagDDEProblem(f, h, [1.0], lags, (0.0, 10.0); iip=DiffEqBase.isinplace(f, 4))
dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test maximum(sol.errors[:l2]) < 1e-4

f = function (t,u,h,du)
    du[1] = - h(t-1)[1]
end
h = (t) -> [0.0]
function (p::typeof(f))(::Type{Val{:analytic}}, t, u0)
    t = t-1
    if t<0
        return u0
    elseif t<1
        return u0 - t*u0
    elseif t<2
        return 0.5*(3u0 -4*t*u0 + t^2 * u0)
    elseif t<3
        return 1/6*(17u0 -24t*u0 + 9t^2*u0 -t^3 * u0)
    elseif t<4
        return 1/24*(149u0 - 204*t*u0 + 90t^2*u0 - 16t^3*u0 + t^4*u0)
    elseif t<5
        return 1/120*(1769u0 - 2300t*u0 + 1090t^2*u0 - 240t^3*u0 + 25t^4*u0 - t^5*u0)
    elseif t<6
        return 1/720*(26239u0 - 32550t*u0 + 15915t^2*u0 - 3940t^3*u0 + 525t^4*u0 - 36t^5*u0 + t^6*u0)
    elseif t<7
        return (463609u0 - 554442t*u0 + 274701t^2*u0 - 72940t^3*u0 + 11235t^4*u0 - 1008t^5*u0 + 49t^6*u0 - t^7*u0)/5040
    elseif t<8
        return (9473673u0 - 11023880t*u0 + 5491780t^2*u0 - 1524712t^3*u0 + 257950t^4*u0 - 27272t^5*u0 + 1764t^6*u0 - 64t^7*u0 + t^8*u0)/40320
    elseif t<9
        return (219480785u0 - 250209864t*u0 + 124923492t^2*u0 - 35742504t^3*u0 +6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
    elseif t<=10
        return (219480785*u0 - 250209864*t*u0 + 124923492*t^2*u0 - 35742504t^3*u0 + 6450318t^4*u0 - 761544t^5*u0 + 58884t^6*u0 - 2880t^7*u0 + 81t^8*u0 - t^9*u0)/362880
    else
        error("This analytical solution is only valid on [-infty,11]")
    end
end
prob = ConstantLagDDEProblem(f, h, [1.0], lags, (0.0, 10.0); iip=DiffEqBase.isinplace(f, 4))
dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test maximum(sol.errors[:l∞]) < 1e-4

lags = [1//3, 1//5]
f = function (t,u,h)
    - h(t-1/3) - h(t-1/5)
end
h = (t) -> 0.0

function (p::typeof(f))(::Type{Val{:analytic}}, t, u0)
    if t<1/5
        return u0
    elseif t<1/3
        return (1/5)*(6u0 - 5t*u0)
    elseif t<2/5
        return (1/15)*(23u0 - 30t*u0)
    elseif t<8/15
        return (1/150)*(242u0 - 360t*u0 + 75t^2*u0)
    elseif t<3/5
        return (1/450)*(854u0 - 1560t*u0 + 675t^2*u0)
    elseif t<2/3
        return (4351u0 - 8205t*u0 + 4050t^2*u0 - 375t^3*u0)/2250
    elseif t<11/15
        return (1/750)*(1617u0 - 3235t*u0 + 1725t^2*u0 - 125t^3*u0)
    elseif t<4/5
        return (7942u0 - 17280t*u0 + 11475t^2*u0 - 2250t^3*u0)/3375
    elseif t<13/15
        return (319984u0 - 702720t*u0 + 480600t^2*u0 - 108000t^3u0 + 5625*t^4u0)/135000
    elseif t<14/15
        return (40436u0 - 94980t*u0 + 72900t^2*u0 - 19500t^3*u0 + 625t^4*u0)/15000
    elseif t<=1
        return (685796u0 - 1670388t*u0 + 1392660t^2*u0 - 467100t^3*u0 + 50625t^4*u0)/243000
    else
        error("This analytical solution is only valid on [-infty,1]")
    end
end

prob = ConstantLagDDEProblem(f, h, 1.0, lags, (0.0, 1.0); iip=DiffEqBase.isinplace(f, 4))
alg = MethodOfSteps(BS3(); constrained=true)

dde_int = init(prob, alg; dt=0.1)
sol = solve!(dde_int)

@test maximum(sol.errors[:final]) < 1e-5
