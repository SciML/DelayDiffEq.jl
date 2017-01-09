# DelayDiffEq.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/DelayDiffEq.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/DelayDiffEq.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ool429apgp28p71x?svg=true)](https://ci.appveyor.com/project/ChrisRackauckas/delaydiffeq-jl)
[![Coverage Status](https://coveralls.io/repos/JuliaDiffEq/DelayDiffEq.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaDiffEq/DelayDiffEq.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaDiffEq/DelayDiffEq.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaDiffEq/DelayDiffEq.jl?branch=master)

OrdinaryDiffEq.jl is a component package in the DifferentialEquations ecosystem. It holds the
delay differential equation solvers and utilities. It is built on top of OrdinaryDiffEq
to extend those solvers for differential delay equations. While completely independent
and usable on its own, users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).
