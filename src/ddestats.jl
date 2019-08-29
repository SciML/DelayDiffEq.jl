mutable struct DDEStats
  nf::Int
  nf2::Int
  nw::Int
  nsolve::Int
  njacs::Int
  nnonliniter::Int
  nnonlinconvfail::Int
  nfpiter::Int
  nfpconvfail::Int
  ncondition::Int
  naccept::Int
  nreject::Int
  maxeig::Float64
end

DDEStats(x::Int = -1) = DDEStats(x, x, x, x, x, x, x, x, x, x, x, x, 0.0)

function Base.show(io::IO, s::DDEStats)
  println(io, summary(s))
  @printf io "%-60s %-d\n" "Number of function 1 evaluations:" s.nf
  @printf io "%-60s %-d\n" "Number of function 2 evaluations:" s.nf2
  @printf io "%-60s %-d\n" "Number of W matrix evaluations:" s.nw
  @printf io "%-60s %-d\n" "Number of linear solves:" s.nsolve
  @printf io "%-60s %-d\n" "Number of Jacobians created:" s.njacs
  @printf io "%-60s %-d\n" "Number of nonlinear solver iterations:" s.nnonliniter
  @printf io "%-60s %-d\n" "Number of nonlinear solver convergence failures:" s.nnonlinconvfail
  @printf io "%-60s %-d\n" "Number of fixed-point solver iterations:" s.nfpiter
  @printf io "%-60s %-d\n" "Number of fixed-point solver convergence failures:" s.nfpconvfail
  @printf io "%-60s %-d\n" "Number of rootfind condition calls:" s.ncondition
  @printf io "%-60s %-d\n" "Number of accepted steps:"  s.naccept
  @printf io "%-60s %-d" "Number of rejected steps:" s.nreject
  iszero(s.maxeig) || @printf io "\n%-60s %-d" "Maximum eigenvalue recorded:" s.maxeig
end