# Defines output directory

# Loads packages
using Plots
using DelimitedFiles
using Statistics
using Images, ImageView
using DifferentialEquations
using DataFrames
using CSV
using Printf

# Include model
include("model_STD_CC_rel_checkFD_SJO.jl")

# Simulation parameters



const T = 1000000
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Tlow = convert(Int64, 1 / dt)
const t = range(dt, T, length = Tdt)

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
const IstepE2 = -0.1
const tstepEinit = T
const tstepEfinal = T

##Cellule I
const IappI = 0.
const IstepI1 = 5.
const IstepI2 = -1.4
const IstepI3 = 5.
const IstepI4 = 5.
const tstepIinit1 = 0
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
const gEC_const = 0.01
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

## Param Deperrois (STD)

const tau_Ca = 38.3492083   # 22.27212 #[ms]
const C_Pre  = 3.99132241   # 0.84410
const C_Post = 1.12940834    # 1.62138
const D_pre  = 9.23545841   # 9.53709 #[ms]

const tau_w    = 299.8778e3 # 520.76129e3 #[s>ms]
const gamma_p  = 564.392975 # 597.08922
const gamma_d  = 111.320539 # 137.7586
const n_lin = 1.
const theta_p = 1.63069609  # 2.009289
const theta_d = 1.
#const eta_lin = ((n_lin*(C_Pre+C_Post)-C_Post)/C_Pre )-1
const U = 0.3838
const tau_rec = 148.9192

## Parameters for Shouval Rule

expm="CC_rel"

a1 = 1.1 #1.
a2 = 1.55
m1 = 0.25
b1 = 40
b2 = 10#40
m2 = 0.5

P1=0.1*4e4
P2=P1*1e-6
P3=2.4
P4=1

#=
a0 = 0.5;
a1 = 1.1;
a2 = 1.47;
b1 = 40;
b2 = 40;
m2 = gamma_p/(gamma_d + gamma_p)


P1=0.1*2.5e4;
P2=P1*0.125e-4;
P3=4.1;
P4=80;
=#
w_init=0.5

const gEC_const = 0.01
wMAX = 1.


##

const count_fMIN = 60
const count_f = 60

frequencies = [1000/0.1,1000/1, 1000/5, 1000/10,1000/15, 1000/20, 1000/25, 1000/30, 1000/35, 1000/40,1000/45, 1000/50]
#frequencies = [1000/0.1,1000/10, 1000/20,  1000/30,  1000/40, 1000/50]

delay_pre = 10.0#-D_pre

wfinal_pos = zeros(length(frequencies))
wfinal_neg = zeros(length(frequencies))
# RUN MODEL
for idx_loop=1:length(frequencies)
    freq=frequencies[idx_loop]
    @time (w, count_pre, count_post, idx_stop) = simulateTOY_ncells(
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
    wfinal_pos[idx_loop] = w[idx_stop]
end

delay_pre = -10.0#-D_pre
# RUN MODEL
for idx_loop=1:length(frequencies)
    freq=frequencies[idx_loop]
    @time (w, count_pre, count_post, idx_stop) = simulateTOY_ncells(
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
    wfinal_neg[idx_loop] = w[idx_stop]
end


wdelta_pos = wfinal_pos/w_init
wdelta_neg = wfinal_neg/w_init

writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF%d_%s_wpos.dat", expm), count_f, wdelta_pos, header=false)
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF%d_%s_wneg.dat", expm), count_f, wdelta_neg, header=false)

plot( [0.1; 1; 5; 10; 15; 20; 25; 30; 35; 40; 45;  50],wdelta_pos, color = "red", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))
plot!([0.1; 1; 5; 10; 15; 20; 25; 30; 35; 40; 45;  50],wdelta_neg, color = "blue", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))

#plot( [0.1; 10; 20; 30; 40; 50],wdelta_pos, color = "red", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))
#plot!([0.1; 10; 20; 30; 40; 50],wdelta_neg, color = "blue", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))
