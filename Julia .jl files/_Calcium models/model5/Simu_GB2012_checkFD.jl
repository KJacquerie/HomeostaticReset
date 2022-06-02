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
include("model_GB2012_checkFD.jl")

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

const gIEGABAB = 1.8 /nIcells
const gIIGABAB = 0.0 /nIcells
const gICGABAB = 1.8 /nIcells


## Param GRAUPNER 2012

#=
const expm="DP_pink"
const strength = "weak"
if(strength=="strong")
    const gEC_const = 0.4
else
    const gEC_const = 0.01
end
const tau_Ca = 20. #[ms]
const C_Pre  = 1.
const C_Post = 2.
const D_pre  = 13.8 #[ms]

const tau_w    = 150e3 #[s>ms]
const gamma_p  = 241.356
const gamma_d  = 150.
const n_lin = 1.
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5
=#


const expm="DP"
const gEC_const = 0.01
const tau_Ca = 20. #[ms]
const C_Pre  = 1.
const C_Post = 2.
const D_pre  = 13.7 #[ms]

const tau_w    = 150e3 #[s>ms]
const gamma_p  = 321.808
const gamma_d  = 200.
const n_lin = 1.
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5

const w_init = 0.5


#=

const expm="cortex"
const gEC_const = 0.01

const tau_Ca = 22.6936 #[ms]
const C_Pre  = 0.5617539
const C_Post = 1.23964
const D_pre  = 4.6098 #[ms]

const tau_w    = 346.3615e3 #[s>ms]
const gamma_p  = 725.085
const gamma_d  = 331.909
#const n_lin = 1.
const theta_p = 1.3
const theta_d = 1.
const wfix = 0.5


w_init = 0.5
=#
const gEC_const = 0.01





##

const count_fMIN = 75
const count_f = 75

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

writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/checkF%d_%s_wpos.dat",count_fMIN, expm), wdelta_pos, header=false
writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/checkF%d_%s_wneg.dat",count_fMIN, expm), wdelta_neg, header=false)

plot( [0.1; 1; 5; 10; 15; 20; 25; 30; 35; 40; 45;  50],wdelta_pos, color = "red", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))
plot!([0.1; 1; 5; 10; 15; 20; 25; 30; 35; 40; 45;  50],wdelta_neg, color = "blue", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))

#plot( [0.1; 10; 20; 30; 40; 50],wdelta_pos, color = "red", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))
#plot!([0.1; 10; 20; 30; 40; 50],wdelta_neg, color = "blue", xlabel = "rho (Hz)", ylabel = " delta w (-)", legend= (frameon=false))



#=

VplotI = plot(
    t[Tlow:Tdt],
    Vconnect[Tlow:Tdt, nEcells+1:ncells-nCcells],
    color = "red",
    ylabel="I [mV]",
    legend= (frameon=false)
)

Vplot1 = plot(
    t[Tlow:Tdt],
    Vconnect[Tlow:Tdt, 1:nEcells],
    ylabel="E [mV]",
    legend= (frameon=false),
    color = "blue")

Vplot2 = plot(
    t[Tlow:Tdt],
    Vconnect[Tlow:Tdt, nEcells+nIcells+1:ncells],
    color = "grey",
    ylabel="C [mV]",
    legend= (frameon=false)
)

cpre_plot = plot(
    t[Tlow:Tdt],
    cpre[Tlow:Tdt],
    color = "blue",
    ylabel="Cpre [-]",
    legend= (frameon=false)
)
cpost_plot = plot(
    t[Tlow:Tdt],
    cpost[Tlow:Tdt],
    color = "gray",
    ylabel="Cpost [-]",
    legend= (frameon=false)
)

cvar_plot = plot(
    t[Tlow:Tdt],
    c_var[Tlow:Tdt],
    color = "black",
    ylabel="C [-]",
    legend= (frameon=false)
)

wplot = plot(
    t[Tlow:Tdt],
    w[Tlow:Tdt],
    color = "black",
    ylabel="w [-]",
    legend= (frameon=false)
)


plot(Vplot1, Vplot2, VplotI, layout = (3, 1))

plot(Vplot1, Vplot2, VplotI,wplot, layout = (4, 1))

plot(cpre_plot, cpost_plot, cvar_plot,wplot, layout=(4,1))
=#
#plot(t[convert(Int64,9000/dt):convert(Int64,10000/dt)], w[convert(Int64,9000/dt):convert(Int64,10000/dt)])
#=
T1=6080
T2=6500
plot(t[convert(Int64,T1/dt):convert(Int64,T2/dt)], Vconnect[convert(Int64,T1/dt):convert(Int64,T2/dt),3])
plot!(t[convert(Int64,T1/dt):convert(Int64,T2/dt)], Vconnect[convert(Int64,T1/dt):convert(Int64,T2/dt),1])
=#
