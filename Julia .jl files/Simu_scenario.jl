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
using Distributions

# Include model

include("model_scenario.jl")

# Simulation parameters

const T = 160000
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


const A2_m = 0.00826477#7e-3
const A3_m = 0.#2.3e-4 #0
const A2_p = 0.#5e-10 #0
const A3_p = 0.0165746#6.5e-3

const tau_m = 33.7
const tau_p = 16.8
const tau_x = 101#101#101
const tau_y = 56.3824 #125


const wMax = 1.
const wMin = 0.

const IappE = 0.
# Step fréquentiel
const IstepE = 35.
const duration = 3.
#const delay_pre = -0.
#Step manuel
const IstepE2 = -0.
const tstepEinit = T
const tstepEfinal = T

const IappI = 3.
const IstepI1 = 0.
const IstepI2 = -IappI-1.2#-IappI-1.2
const IstepI3 = 0.
const IstepI4 = -IappI-1.3#-IappI-1.3
const tstepIinit1 = 0
const tstepIinit2 = Int64(T/4)
const tstepIinit3 = 2*Int64(T/4)
const tstepIinit4 = 3*Int64(T/4)
const tstepIfinal = copy(T)

const IappC = 0.
#Step fréquentiel
const IstepC = 35.
#const duration = 20.0
const IstepC2 = -0.
const tStartC = T
const tStopC = T

# connectivity
#const gEEAMPA = 0.0/nEcells

const gEC_const = 0.01
const gEI = 0.0 / nEcells
const gEE = 0.0 / nEcells

const gIEGABAA = 2. / nIcells
const gICGABAA = 2. / nIcells
const gIIGABAA = 0.0 / nIcells

const gCE = 0.0 / nCcells
const gCI = 0.0 / nCcells
const gCC = 0.0 / nCcells

const gIEGABAB = 1. /nIcells
const gIIGABAB = 0.0 /nIcells
const gICGABAB = 1. /nIcells
w_init = 0.5

# freq_ = 50.
# freq_wake = 1/freq_/0.001
delay_pre = 5.

##
#
# freq_ = 30.
# println(freq_)
# freq_wake = round(1/freq_/0.001)
#
# @time ( Vconnect, w ) = debug_scenario_conso(
#             freq_wake,
#             delay_pre,
#             w_init)
#
# VplotI = plot(
#     t[Tlow:Tdt],
#     Vconnect[Tlow:Tdt, nEcells+1:ncells-nCcells],
#     color = "red",
#     ylabel="I [mV]",
#     legend= (frameon=false)
# )
#
# Vplot1 = plot(t[Tlow:Tdt], Vconnect[Tlow:Tdt, 1:nEcells],ylabel="E [mV]", legend= (frameon=false)
# , color = "blue")
#
# Vplot2 = plot(
#     t[Tlow:Tdt],
#     Vconnect[Tlow:Tdt, nEcells+nIcells+1:ncells],
#     color = "grey",
#     ylabel="C [mV]",
#     legend= (frameon=false)
# )
#
# Vplot3 = plot(t[Tlow:Tdt], w[Tlow:Tdt], ylabel="w(t) [-]", color = "green", legend= (frameon=false))
#
# plot(Vplot1, Vplot2, VplotI, Vplot3, layout = (4, 1))
#

##
nb_circuits= 6
freq_vec =[ 5. , 10., 20. ]
percentage_conso = 0.5

cd("C:/Users/carol.000/Documents/PhD/Data/scenario")

for i = 4:nb_circuits

    if i<= nb_circuits*percentage_conso
        freq_ = freq_vec[i]
        #freq_ = rand([30., 40., 50.])
        println(freq_)
        freq_wake = round(1/freq_/0.001)
        @time ( w , V) = simulate_scenario_conso(
                    freq_wake,
                    delay_pre,
                    w_init)
    else
        freq_ = freq_vec[i-3]
        #freq_ = rand([30., 40., 50.])
        #freq_ = rand(30:50)
        freq_wake = round(1/freq_/0.001)
        @time ( w, V) = simulate_scenario_noisy(
                    freq_wake,
                    delay_pre,
                    w_init)
    end

    writedlm(@sprintf("w_neuron%d.dat", i), w, header = false)
    writedlm(@sprintf("V_E%d.dat", i), V[:,1], header = false)
    writedlm(@sprintf("V_C%d.dat", i), V[:,3], header = false)

end


## Kernel
# deltas = range(-50, 50, length= 10)
#
# wSTDP = zeros(length(deltas))
#
# freq = 100.

# w_init = 0.3
#
# for i = 1:length(deltas)
#
#         delay_pre = deltas[i]
#         @time ( w_ ) = simulate_regular(
#             ncells,
#             nEcells,
#             nIcells,
#             IappE,
#             IappI,
#             tstepEinit,
#             tstepEfinal,
#             IstepE,
#             tstepIinit1,
#             tstepIinit2,
#             tstepIinit3,
#             tstepIinit4,
#             tstepIfinal,
#             IstepI1,
#             IstepI2,
#             IstepI3,
#             IstepI4,
#             gEC_const,
#             freq,
#             delay_pre,
#             w_init
#         )
#         wSTDP[i] = w_/w_init#( w_ - w_init) *100 / w_init
#
# end
# #
# plot(deltas, wSTDP ,color = "blue", xlims = (-50, 50), legend = false)

## RUN MODEL: time course

# delay_pre = 10.
# @time (Vconnect, w, x_, y_) = simulate_timecourse(
#     ncells,
#     nEcells,
#     nIcells,
#     IappE,
#     IappI,
#     tstepEinit,
#     tstepEfinal,
#     IstepE,
#     tstepIinit1,
#     tstepIinit2,
#     tstepIinit3,
#     tstepIinit4,
#     tstepIfinal,
#     IstepI1,
#     IstepI2,
#     IstepI3,
#     IstepI4,
#     gEC_const,
#     freq,
#     delay_pre,
#     w_init
# )
# # #println(stop/dt)
# # #println(w_[Int64(round(stop/dt))]-w_[1])
# #
# # VplotI = plot(
# #     t[Tlow:Tdt],
# #     Vconnect[Tlow:Tdt, nEcells+1:ncells-nCcells],
# #     color = "red",
# #     ylabel="I [mV]",
# #     legend= (frameon=false)
# # )
# #
# # Vplot1 = plot(t[Tlow:Tdt], Vconnect[Tlow:Tdt, 1:nEcells],ylabel="E [mV]", legend= (frameon=false)
# # , color = "blue")
# #
# # Vplot2 = plot(
# #     t[Tlow:Tdt],
# #     Vconnect[Tlow:Tdt, nEcells+nIcells+1:ncells],
# #     color = "grey",
# #     ylabel="C [mV]",
# #     legend= (frameon=false)
# # )
# #
# # Vplot3 = plot(t[Tlow:Tdt], w[Tlow:Tdt], ylabel="w(t) [-]", color = "green", legend= (frameon=false))
# # Vplotx = plot(t[Tlow:Tdt], x_[Tlow:Tdt], ylabel="x(t) [-]", color = "green", legend= (frameon=false))
# # Vploty = plot(t[Tlow:Tdt], y_[Tlow:Tdt], color = "green",ylabel="y(t) [-]", legend= (frameon=false))
# #
# #
# # plot(Vplot1, Vplot2, Vplot3, Vplotx, Vploty, layout = (5, 1))
#
# # cd(@sprintf("C:/Users/carol.000/Documents/PhD/Data/temporel/w_%s", strength))
# # writedlm("tBurst.dat", t, header = false)
# # writedlm("wBurst.dat", w, header = false)
# # writedlm("VBurst.dat", Vconnect, header = false)
# # writedlm("xBurst.dat", x_, header = false)
# # writedlm("yBurst.dat", y_, header = false)
