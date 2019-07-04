using DelayDiffEq, Test
using DiffEqProblemLibrary.DDEProblemLibrary

DDEProblemLibrary.importddeproblems()

using DiffEqProblemLibrary.DDEProblemLibrary:
  # DDE problems with 1 constant delay
  prob_dde_constant_1delay_ip, prob_dde_constant_1delay_oop, prob_dde_constant_1delay_scalar,
  prob_dde_constant_1delay_long_ip, prob_dde_constant_1delay_long_oop, prob_dde_constant_1delay_long_scalar,
  # DDE problems with 2 constant delays
  prob_dde_constant_2delays_ip, prob_dde_constant_2delays_oop, prob_dde_constant_2delays_scalar,
  prob_dde_constant_2delays_long_ip, prob_dde_constant_2delays_long_oop, prob_dde_constant_2delays_long_scalar,
  # Mackey-Glass equation (model of blood production) with 1 constant delay
  prob_dde_DDETST_A1
