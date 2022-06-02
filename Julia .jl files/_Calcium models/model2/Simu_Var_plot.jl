# Defines output directory

# Loads packages
using Plots
using DelimitedFiles
using Statistics
using Images, ImageView
using DifferentialEquations
using DataFrames
using Printf
using CSV

# Include model
include("model_Control_plot.jl")

# Simulation parameters

const T = 6000
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Tlow = convert(Int64, 1 / dt)
const t = range(dt, T, length = Tdt)

const crit_burst=3.

# Model parameters (global)
const C = 1
const VNa = 50
const VK = -85
const VCa = 120
const Vl = -55
const VH = -20
const Kd = 170


# Model parameters (mean) - Ecells
const gl_E = 0.055
const gNa_E = 170
const gKd_E = 40
const k1_E = 1.e-1
const k2_E = 0.1e-1
const gH_E = 0.01
const gKCa_E = 4#/10
const gCaT_E = 0.55

# Model parameters (mean) - Icells
const gl_I = 0.055
const gNa_I = 170
const gKd_I = 40
const k1_I = 1.e-1
const k2_I = 0.1e-1
const gH_I = 0.01
const gKCa_I = 4#4#/10
const gCaT_I = 0.55

# Model parameters (mean) - Ccells
const gl_C = 0.055
const gNa_C = 170
const gKd_C = 40
const k1_C = 1.e-1
const k2_C = 0.1e-1
const gH_C = 0.01

const idx_g=2
if idx_g == 1
    const gKCa_C = 4.#*0.75
    const gCaT_C = 0.55#*1.5
    const C_Post = 1.62138
else
    const gKCa_C = 4*0.75
    const gCaT_C = 0.55*1.5
    const C_Post = 1.62138*1.5
end





# Simulations
const nEcells = 1 # Number of excitatory cells
const nIcells = 1 # Number of inhibitory cells
const nCcells = 1 # Number of cortical cells
const ncells = nEcells + nIcells + nCcells

const gNavec_E = gNa_E * ones(ncells)
const gKdvec_E = gKd_E * ones(ncells)
const gHvec_E = gH_E * ones(ncells)
const gCaTvec_E = gCaT_E * ones(ncells)
const gKCavec_E = gKCa_E * ones(ncells)
const glvec_E = gl_E * ones(ncells)
const k1vec_E = k1_E * ones(ncells)
const k2vec_E = k2_E * ones(ncells)

const gNavec_I = gNa_I * ones(ncells)
const gKdvec_I = gKd_I * ones(ncells)
const gHvec_I = gH_I * ones(ncells)
const gCaTvec_I = gCaT_I * ones(ncells)
const gKCavec_I = gKCa_I * ones(ncells)
const glvec_I = gl_I * ones(ncells)
const k1vec_I = k1_I * ones(ncells)
const k2vec_I = k2_I * ones(ncells)

const gNavec_C = gNa_C * ones(ncells)
const gKdvec_C = gKd_C * ones(ncells)
const gHvec_C = gH_C * ones(ncells)
const gCaTvec_C = gCaT_C * ones(ncells)
const gKCavec_C = gKCa_C * ones(ncells)
const glvec_C = gl_C * ones(ncells)
const k1vec_C = k1_C * ones(ncells)
const k2vec_C = k2_C * ones(ncells)


##Cellule E
const IappE = 0.
# Step fréquentiel
const IstepE = 50.
const duration = 3.
#const delay_pre = -0.
#Step manue
const IstepE2 = 0.#-0.1
const tstepEinit = 1000
const tstepEfinal = T

##Cellule I
const IappI = 3.#0.
const IstepI1 = -1.2-IappI#5.
#=
if (IstepI1 ==-4.2)
    const idx_I=2
else
    const idx_I=1
end
=#
const IstepI2 = 0.#-1.4
const IstepI3 = 0.#5.
const IstepI4 = 0.#5.
const tstepIinit1 = tstepEinit
const tstepIinit2 = T
const tstepIinit3 = T
const tstepIinit4 = T
const tstepIfinal = copy(T)


##Cellule C
const IappC = IappE
#Step fréquentiel
const IstepC = 50.
#const duration = 20.0
const IstepC2 = IstepE2
const tStartC = tstepEinit
const tStopC = T


# connectivity

const gEI = 0.0 / nEcells
const gEE = 0.0 / nEcells

const gIEGABAA = 2. / nIcells
const gICGABAA = 2. / nIcells
const gIIGABAA = 0.0 / nIcells

const gCE = 0.0 / nCcells
const gCI = 0.0 / nCcells
const gCC = 0.0 / nCcells

const gIEGABAB = 1.5/nIcells
const gIIGABAB = 0.0 /nIcells
const gICGABAB = 1.5 /nIcells



## Param GRAUPNER 2016


const expm = "Control"
const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410

const D_pre  = 9.53709 #[ms]

const tau_w    = 520.76129e3 #[s>ms]
const gamma_p  = 597.08922
const gamma_d  = 137.7586
const n_lin = 1.
const theta_p = 2.009289
const theta_d = 1.0


Omega_p = gamma_p/(gamma_p+gamma_d)
Omega_d = 0.
Omega_0 = 0.

tauw_p = tau_w /(gamma_p + gamma_d)
tauw_d = tau_w / gamma_d
tauw_0 = 0.

const w_init = 0.5
const strength = "weak"
if(strength=="strong")
    const gEC_const = 0.4
else
    const gEC_const = 0.01
end




##

freq_intensity=10.
freq = 1000/freq_intensity
delay_pre = -10.0#-D_pre#-10.0-D_pre
# RUN MODEL
@time (Vconnect, cpre, cpost, c_var, spk_pre, spk_post, w, CaT) = simulateTOY_ncells(
    ncells,
    nEcells,
    nIcells,
    IappE,
    IappI,
    tstepEinit,
    tstepEfinal,
    IstepE,
    tstepIinit1,
    tstepIinit2,
    tstepIinit3,
    tstepIinit4,
    tstepIfinal,
    IstepI1,
    IstepI2,
    IstepI3,
    IstepI4,
    gEC_const,
    freq,
    delay_pre
)



T1 = 1#250000
T2 = Int64(T/dt)#300000



VplotI = plot(
    t[T1:T2],
    Vconnect[T1:T2,2],
    color = "red",
    ylabel="I [mV]",
    legend= (frameon=false)
)

Vplot1 = plot(
    t[T1:T2],
    Vconnect[T1:T2,1],
    color = "blue",
    ylabel="E [mV]",
    legend= (frameon=false)
)

Vplot2 = plot(
    t[T1:T2],
    Vconnect[T1:T2,3],
    color = "grey",
    ylabel="C [mV]",
    legend= (frameon=false)
)



cpre_plot = plot(
    t[T1:T2],
    cpre[T1:T2],
    color = "blue",
    ylabel="Cpre [-]",
    legend= (frameon=false)
)
cpost_plot = plot(
    t[T1:T2],
    cpost[T1:T2],
    color = "gray",
    ylabel="Cpost [-]",
    legend= (frameon=false)
)

cvar_plot = plot(
    t[T1:T2],
    c_var[T1:T2],
    color = "black",
    ylabel="C [-]",
    legend= (frameon=false)
)

wplot = plot(
    t[T1:T2],
    w[T1:T2],
    color = "black",
    ylabel="w [-]",
    legend= (frameon=false)
)



writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/Vconnect_g%d.dat",idx_g), Vconnect, header = false)
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/w_g%d.dat",idx_g), w, header = false)
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/cpre_g%d.dat",idx_g), cpre, header = false)
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/cpost_g%d.dat",idx_g), cpost, header = false)
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/ca_g%d.dat",idx_g), c_var, header = false)


plot(Vplot1, Vplot2, VplotI,wplot, layout = (4, 1))
