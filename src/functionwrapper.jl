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
  if f.jac === nothing
    jac = nothing
  else
    if isinplace(f)
      jac = let f_jac = f.jac, h = h
        (J, u, p, t) -> f_jac(J, u, h, p, t)
      end
    else
      jac = let f_jac = f.jac, h = h
        (u, p, t) -> f_jac(u, h, p, t)
      end
    end
  end

  ODEFunctionWrapper{isinplace(f),typeof(f.f),typeof(h),typeof(f.mass_matrix),
                     typeof(f.analytic),typeof(f.tgrad),typeof(jac),
                     typeof(f.jac_prototype),typeof(f.Wfact),typeof(f.Wfact_t),
                     typeof(f.paramjac),typeof(f.syms),typeof(f.colorvec)}(
                       f.f, h, f.mass_matrix, f.analytic, f.tgrad, jac,
                       f.jac_prototype, f.Wfact, f.Wfact_t, f.paramjac, f.syms,
                       f.colorvec)
end

(f::ODEFunctionWrapper{true})(du, u, p, t) = f.f(du, u, f.h, p, t)
(f::ODEFunctionWrapper{false})(u, p, t) = f.f(u, f.h, p, t)
