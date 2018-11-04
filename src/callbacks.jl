# update integrator when u is modified by callbacks
function handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities, Discontinuity(integrator.t, 0))
end

"""
    reeval_internals_due_to_modification!(integrator::DDEIntegrator)

Recalculate interpolation data and update ODE integrator after changes by callbacks.
"""
function DiffEqBase.reeval_internals_due_to_modification!(integrator::DDEIntegrator,
            x::Type{Val{not_initialization}} = Val{true}) where not_initialization

    # update interpolation data of DDE integrator using old interpolation data
    # of ODE integrator in evaluation of history function that was calculated in
    # `perform_step!`
    if not_initialization
        addsteps!(integrator, integrator.f, true, true, true)
        # copy interpolation data to ODE integrator
        recursivecopy!(integrator.integrator.k, integrator.k)
    end

    # move ODE integrator to new time interval of DDE integrator
    integrator.integrator.t = integrator.t
    integrator.integrator.dt = integrator.dt
    if isinplace(integrator.sol.prob)
        recursivecopy!(integrator.integrator.u, integrator.u)
    else
        integrator.integrator.u = integrator.u
    end

    integrator.u_modified = false
end

# track propagated discontinuities for dependent delays

struct DiscontinuityCallback{F,D<:Discontinuity,A,R,I} <: AbstractContinuousCallback
    lags::F
    discontinuities::Vector{D}
    interp_points::Int
    abstol::A
    reltol::R
    initialize::I
    idxs::Nothing
end

"""
    DiscontinuityCallback(lags, discontinuities::Vector{<:Discontinuity};
                          [interp_points::Int=10, abstol=1e-12, reltol=0])

Callback that tracks `discontinuities` that are propagated by dependent `lags` of the form
`(u,p,t) -> lag(u,p,t)`.

Hereby a number `interp_points` of interpolation points are used to first check for
different signs of functions ``f(t) = T + lag(u(t),p,t) - t``, where ``T`` is time point of a
previous discontinuity and ``t`` is contained in the current time interval. This shows that
the current time interval contains propagated discontinuities of which the exact time point
is then determined by a root finding algorithm. The sign at the lower bound of the time
interval, i.e. at `tprev`, is set to 0 with absolute tolerance `abstol` and relative
tolerance `reltol`.
"""
function DiscontinuityCallback(lags, discontinuities::Vector{<:Discontinuity};
                               interp_points::Int=10, abstol=1e-12, reltol=0)
    DiscontinuityCallback(lags, discontinuities, interp_points, abstol,
                          reltol, initialize!, nothing)
end

# do not initialize discontinuity callback
initialize!(c::DiscontinuityCallback, u, t, integrator::DEIntegrator) = (integrator.u_modified=false)

# find time of first discontinuity in the current time interval (if existent)
function DiffEqBase.find_callback_time(integrator::DDEIntegrator, callback::DiscontinuityCallback, counter)
    # initialize time and order of first discontinuity in the current time interval
    tmin = zero(integrator.t)
    order = 0
    contains_discontinuity = false

    # update interpolation data and calculate interpolation points
    if callback.interp_points != 0
      addsteps!(integrator)
    end
    Θs = range(0, stop=one(integrator.t), length=callback.interp_points)

    for lag in callback.lags
        # define function f that calculates T + lag(u(t),p,t) - t, where t = tprev + Θ*dt is
        # a time point in the time interval of the current step and T ≤ tprev is the time
        # point of a discontinuity
        # hence roots t of this function for fixed T are propagated discontinuities
        function f(Θ, T)
            if typeof(integrator.cache) <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
                # use temporary array of mutable caches in calculation of interpolants
                tmp = integrator.cache.tmp
                OrdinaryDiffEq.ode_interpolant!(tmp, Θ, integrator, nothing, Val{0})
            else
                tmp = OrdinaryDiffEq.ode_interpolant(Θ, integrator, nothing, Val{0})
            end
            t = integrator.tprev + Θ*integrator.dt
            T + lag(tmp,integrator.p,t) - t
        end

        for d in callback.discontinuities
            # fix time of discontinuity
            T = d.t
            g(Θ) = f(Θ, T)

            # use start and end point of last time interval to check for discontinuities
            previous_condition = T + lag(integrator.uprev, integrator.p, integrator.tprev) -
                integrator.tprev
            if isapprox(previous_condition, 0, rtol=callback.reltol, atol=callback.abstol)
                prev_sign = 0
            else
                prev_sign = cmp(previous_condition, 0)
            end
            next_sign = cmp(T + lag(integrator.u, integrator.p, integrator.t) - integrator.t, 0)

            # find time of discontinuity if one exists
            t,contains_discontinuity = find_discontinuity_time(integrator, callback, prev_sign, next_sign, Θs, g)

            # update time and order of first discontinuity in the current time interval
            if contains_discontinuity && (integrator.tdir * t < integrator.tdir * tmin || iszero(tmin))
                tmin = t
                order = d.order + 1
            end
        end
    end

    tmin, order, contains_discontinuity
end

"""
    find_discontinuity_time(integrator::DDEIntegrator, callback::DiscontinuityCallback,
                            prev_sign::Int, next_sign::Int, Θs, f)

Find time of propagated discontinuity for a certain dependent delay and previous
discontinuity, which is root of the function `f`, in the current time interval of
`integrator` with interpolation points `Θs`. Hereby `f` shows signs `prev_sign` and
`next_sign` at both ends of the time interval.
"""
function find_discontinuity_time(integrator::DDEIntegrator, callback::DiscontinuityCallback,
                                 prev_sign::Int, next_sign::Int, Θs, f)
    if prev_sign * next_sign < 0
        # if signs are different we already know that a discontinuity exists
        contains_discontinuity = true
        interp_index = callback.interp_points
        prev_sign_index = 1
    elseif callback.interp_points != 0
        # recheck interpolation intervals if no discontinuity found yet
        contains_discontinuity, interp_index, prev_sign_index =
            determine_discontinuity_existence(prev_sign, next_sign, Θs, f)
    end

    if contains_discontinuity
        if interp_index != 0
            top_Θ = Θs[interp_index] # Top at the smallest
            bottom_θ = Θs[prev_sign_index]
        else
            # without interpolation points
            top_Θ = one(integrator.t)
            bottom_θ = typeof(integrator.t)(0)
        end

        if f(top_Θ) == 0
          Θ = top_Θ
        else
          Θ = prevfloat(find_zero(f,(bottom_θ,top_Θ),Roots.AlefeldPotraShi(),atol = callback.abstol/100))
        end
        #Θ = prevfloat(...)
        # prevfloat guerentees that the new time is either 1 floating point
        # numbers just before the event or directly at zero, but not after.
        # If there's a barrier which is never supposed to be crossed,
        # then this will ensure that
        # The item never leaves the domain. Otherwise Roots.jl can return
        # a float which is slightly after, making it out of the domain, causing
        # havoc.
        new_t = integrator.dt*Θ
    else
        # no discontinuity found for given combination of lag and previous discontinuity
        new_t = zero(integrator.t)
    end

    new_t,contains_discontinuity
end

"""
    determine_discontinuity_existence(prev_sign::Int, next_sign::Int, Θs, f)

Determine whether function `f` has a root in the interval [0, 1] by checking signs of `f`
at 0 and 1 (`prev_sign` and `next_sign`, respectively) and at interpolation points in `Θs`.

This corresponds to the existence of a propagated discontinuity.
"""
function determine_discontinuity_existence(prev_sign::Int, next_sign::Int, Θs, f)
    contains_discontinuity = false
    interp_index = 0
    prev_sign_index = 1

    # use interpolation intervals to check for discontinuities
    for i in 2:length(Θs)-1
        new_sign = cmp(f(Θs[i]), 0)

        if prev_sign == 0
            prev_sign = new_sign
            prev_sign_index = i
        elseif prev_sign * new_sign < 0
            contains_discontinuity = true
            interp_index = i
            break
        end
    end

    contains_discontinuity, interp_index, prev_sign_index
end

"""
    apply_callback!(integrator::DDEIntegrator, callback::DiscontinuityCallback, cb_time,
                    order)

Handle discontinuity of order `order` at time `integrator.tprev + cb_time` in the current
time interval of `integrator` that was found by `callback`. Cause the current step to fail,
and add the found discontinuity to both the heap of discontinuities and of time stops.
"""
function DiffEqBase.apply_callback!(integrator::DDEIntegrator, callback::DiscontinuityCallback,
                         cb_time, order)
    # do not accept current step
    integrator.t = integrator.tprev
    integrator.force_stepfail = true
    integrator.accept_step = false

    # add new discontinuity of correct order at time of callback
    d = Discontinuity(integrator.t + cb_time, order)
    push!(integrator.opts.d_discontinuities, d)
    push!(integrator.opts.tstops, d.t)

    # u is not modified and step not saved
    # nevertheless callback set `save_cb` to `true` to avoid undesired saving in
    # `handle_callbacks!`
    false, true
end
