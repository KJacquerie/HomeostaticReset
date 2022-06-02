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
include("model_Control_plotBurst.jl")

# Simulation parameters

const T = 5000
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Tlow = convert(Int64, 500 / dt)
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

const gCaT_WT =0.55
const gCaT_MUT = 0.55

# Model parameters (mean) - Ecells
const gl_E = 0.055
const gNa_E = 170
const gKd_E = 40
const k1_E = 1.e-1
const k2_E = 0.1e-1
const gH_E = 0.01
const gKCa_E = 4#/10
const gCaT_E = 0.55#/5

# Model parameters (mean) - Icells
const gl_I = 0.055
const gNa_I = 170
const gKd_I = 40
const k1_I = 1.e-1
const k2_I = 0.1e-1
const gH_I = 0.01
const gKCa_I = 4#4#/10
const gCaT_I = 0.55#0.55#/5

# Model parameters (mean) - Ccells
const gl_C = 0.055
const gNa_C = 170
const gKd_C = 40
const k1_C = 1.e-1
const k2_C = 0.1e-1
const gH_C = 0.01
const gKCa_C = 4#4#/10
const gCaT_C = 0.55#0.55#/5

# Simulations
const nEcells = 1 # Number of excitatory cells
const nIcells = 1 # Number of inhibitory cells
const nCcells = 1 # Number of cortical cells
const ncells = nEcells + nIcells + nCcells

const gNavec_E = gNa_E * ones(ncells) #- 70*(rand(ncells)-0.5)
const gKdvec_E = gKd_E * ones(ncells) #- 40*(rand(ncells)-0.5)
const gHvec_E = gH_E * ones(ncells) #- 0.002*(rand(ncells)-0.5)
const gCaTvec_E = gCaT_E * ones(ncells) #- 0.35*(rand(ncells)-0.5)
const gKCavec_E = gKCa_E * ones(ncells) #- 2*(rand(ncells)-0.5)
const glvec_E = gl_E * ones(ncells) #- 0.015*(rand(ncells)-0.5)
const k1vec_E = k1_E * ones(ncells) #- 0.05*(rand(ncells))
const k2vec_E = k2_E * ones(ncells) #- 0.005*(rand(ncells))

const gNavec_I = gNa_I * ones(ncells) #- 70*(rand(ncells)-0.5)
const gKdvec_I = gKd_I * ones(ncells) #- 40*(rand(ncells)-0.5)
const gHvec_I = gH_I * ones(ncells) #- 0.002*(rand(ncells)-0.5)
const gCaTvec_I = gCaT_I * ones(ncells) #- 0.35*(rand(ncells)-0.5)
const gKCavec_I = gKCa_I * ones(ncells) #- 2*(rand(ncells)-0.5)
const glvec_I = gl_I * ones(ncells) #- 0.015*(rand(ncells)-0.5)
const k1vec_I = k1_I * ones(ncells) #- 0.05*(rand(ncells))
const k2vec_I = k2_I * ones(ncells) #- 0.005*(rand(ncells))


##Cellule E
const IappE = 0.
# Step fréquentiel
const IstepE = 50.
const duration = 3.
#const delay_pre = -0.
#Step manue
const IstepE2 = 0.#-0.1
const tstepEinit = 100
const tstepEfinal = T

##Cellule I
const IappI = 3.#0.
#const IstepI1 = -1.2-IappI#5.
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

const gIEGABAB = 1.5 /nIcells
const gIIGABAB = 0.0 /nIcells
const gICGABAB = 1.5 /nIcells



## Param GRAUPNER 2016


const expm = "Control"
const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410
const C_Post = 1.62138
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

w_init = 0.5
const strength = "weak"
if(strength=="strong")
    const gEC_const = 0.4
else
    const gEC_const = 0.01
end




##

freq_intensity=10.
freq = 1000/freq_intensity
delay_pre = 10.0#-D_pre#-10.0-D_pre

const IstepI_vec = -1.7-IappI:0.1:-0.9-IappI
idx_Ivec=9
IstepI1 = IstepI_vec[idx_Ivec]

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


const T1 = convert(Int64, 1 / dt)
const T2 = convert(Int64, 5000 / dt)



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


#plot(Vplot1, Vplot2, VplotI, layout = (3, 1))

plot(Vplot1, Vplot2, VplotI,wplot, layout = (4, 1))
plot(cpre_plot, cpost_plot, cvar_plot,wplot, layout=(4,1))
plot(Vplot1, cpre_plot, Vplot2, cpost_plot, cvar_plot, wplot, layout = (6,1))

plot(t[T1:T2], Vconnect[T1:T2,1], color = "blue", ylabel=" [mV]", legend= (frameon=false))
Vplotcomb = plot!(t[T1:T2], Vconnect[T1:T2,3], color = "gray", ylabel=" [mV]", legend= (frameon=false))

plot(t[T1:T2], cpre[T1:T2], color = "blue", ylabel="Ca [-]", legend= (frameon=false))
Cacomb = plot!(t[T1:T2], cpost[T1:T2], color = "gray", ylabel="Ca [-]", legend= (frameon=false))

plot(Vplotcomb, Cacomb, cvar_plot, wplot, layout = (4,1))


writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model2/trace/burst/Vconnect_%d.dat",idx_Ivec), Vconnect, header = false)
