# update integrator when u is modified by callbacks
function OrdinaryDiffEq.handle_callback_modifiers!(integrator::DDEIntegrator)
    integrator.reeval_fsal = true # recalculate fsalfirst after applying step

    # update heap of discontinuities
    # discontinuity is assumed to be of order 0, i.e. solution x is discontinuous
    push!(integrator.opts.d_discontinuities, Discontinuity(integrator.t, 0))
end

# recalculate interpolation data and update the ODE integrator
function DiffEqBase.reeval_internals_due_to_modification!(integrator::DDEIntegrator,
            x::Type{Val{not_initialization}} = Val{true}) where not_initialization
  ode_integrator = integrator.integrator

  if not_initialization
    # update interpolation data of the integrator using the old dense history
    # of the ODE integrator
    DiffEqBase.addsteps!(integrator, integrator.f, true, true, true)

    # copy interpolation data to the ODE integrator
    recursivecopy!(ode_integrator.k, integrator.k)
    end

  # adjust current interval of the ODE integrator
  ode_integrator.t = integrator.t
  ode_integrator.dt = integrator.dt
  if isinplace(integrator.sol.prob)
    recursivecopy!(ode_integrator.u, integrator.u)
  else
    ode_integrator.u = integrator.u
  end

  integrator.u_modified = false
end
