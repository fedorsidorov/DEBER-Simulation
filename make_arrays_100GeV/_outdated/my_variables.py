#%% Variables
e = 4.8e-10
m = 9.11e-28
m_eV = 511e+3
h = 6.63e-27
eV = 1.6e-12
E0 = 20e+3
Na = 6.02e+23

Z_H = 1
u_H = 1.01
rho_H = 8.988e-5
n_H =  rho_H*Na/u_H

Z_C = 6
u_C = 12
rho_C = 2.265
n_C =  rho_C*Na/u_C

Z_O = 8
u_O = 16
rho_O = 1.429e-3
n_O =  rho_O*Na/u_O

Z_Si = 14
u_Si = 28.08
rho_Si = 2.33
n_Si =  rho_Si*Na/u_Si

u_PMMA = 100
rho_PMMA = 1.18
n_PMMA =  rho_PMMA*Na/u_PMMA
n_PMMA_at = n_PMMA*(5 + 2 + 8)

CONC = [n_PMMA_at, n_PMMA_at, n_PMMA_at, n_Si]