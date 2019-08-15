# convenience macro
macro wrap_h(signature)
  Meta.isexpr(signature, :call) ||
    throw(ArgumentError("signature has to be a function call expression"))

  name = signature.args[1]
  args = signature.args[2:end]
  args_wo_h = [arg for arg in args if arg !== :h]
  
  quote
    if f.$name === nothing
      nothing
    else
      if isinplace(f)
        let _f = f.$name, h = h
          ($(args_wo_h...),) -> _f($(args...))
        end
      else
        let _f = f.$name, h = h
          ($(args_wo_h[2:end]...),) -> _f($(args[2:end]...))
        end
      end
    end
  end |> esc
end

struct ODEFunctionWrapper{iip,F,H,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: DiffEqBase.AbstractODEFunction{iip}
  f::F
  h::H
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  syms::S
  colorvec::TCV
end

function ODEFunctionWrapper(f::DDEFunction, h)
  # wrap functions
  jac = @wrap_h jac(J, u, h, p, t)
  Wfact = @wrap_h Wfact(W, u, h, p, dtgamma, t)
  Wfact_t = @wrap_h Wfact_t(W, u, h, p, dtgamma, t)

  ODEFunctionWrapper{isinplace(f),typeof(f.f),typeof(h),typeof(f.mass_matrix),
                     typeof(f.analytic),typeof(f.tgrad),typeof(jac),
                     typeof(f.jac_prototype),typeof(Wfact),typeof(Wfact_t),
                     typeof(f.paramjac),typeof(f.syms),typeof(f.colorvec)}(
                       f.f, h, f.mass_matrix, f.analytic, f.tgrad, jac,
                       f.jac_prototype, Wfact, Wfact_t, f.paramjac, f.syms,
                       f.colorvec)
end

(f::ODEFunctionWrapper{true})(du, u, p, t) = f.f(du, u, f.h, p, t)
(f::ODEFunctionWrapper{false})(u, p, t) = f.f(u, f.h, p, t)
