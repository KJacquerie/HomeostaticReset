# Defines output directory

# Loads packages
using Plots
using DelimitedFiles
using Statistics
using Images, ImageView
using DataFrames
using CSV
using Distributions
# Include model
include("model_SJO_Triplet.jl")


## Simulation parameters
const T = 1000000
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Tlow = convert(Int64, 1 / dt)
const t = range(dt, T, length = Tdt)

## Model parameters (global)
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
const gKCa_C = 4#/10
const gCaT_C = 0.55#/5

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

## Plasticity parameters
# #Minimal model Pfister2006 Cortex (Hard Bounds)
# const A2_m = 7.1e-3#7e-3
# const A3_m = 0.#2.3e-4 #0
# const A2_p = 0.#5e-10 #0
# const A3_p = 6.5e-3#6.5e-3
#
# const tau_m = 33.7
# const tau_p = 16.8
# const tau_x = 101#101#101
# const tau_y = 114.#125
# SB = 0

# #Minimal model Pfister2006 Hippocampe (Hard bounds or soft bounds (SB=1))
# const A2_m = 3.5e-3#7e-3
# const A3_m = 0.#2.3e-4 #0
# const A2_p = 5.3e-3#5e-10 #0
# const A3_p = 8e-3#6.5e-3
#
# const tau_m = 33.7
# const tau_p = 16.8
# const tau_x = 101#101#101
# const tau_y = 40.#125
# SB = 0

#Graupner2016 parameters Cortex (Soft bounds)
# const A2_m = 0.00826477#7e-3
# const A3_m = 0.#2.3e-4 #0
# const A2_p = 0.#5e-10 #0
# const A3_p = 0.0165746#6.5e-3
#
# const tau_m = 33.7
# const tau_p = 16.8
# const tau_x = 101#101#101
# const tau_y = 56.3824 #125
# SB = 1

const wMax = 2.
const wMin = 0.

const IappE = 0.
const IstepE = 50.

# Step fr??quentiel
const duration = 5.
#Step manue
const IstepE2 = 0.
const tstepEinit = T
const tstepEfinal = T

const IappI = 0.
const IstepI1 = 5.
const IstepI2 = -1.
const IstepI3 = 5.
const IstepI4 = -1.
const tstepIinit1 = 0
const tstepIinit2 = T
const tstepIinit3 = T
const tstepIinit4 = T
const tstepIfinal = copy(T)

const IappC = 0.
const IstepC = 50.

const IstepC2 = 0.
const tStartC = T
const tStopC = T

## connectivity
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

const gIEGABAB = 1.5 /nIcells
const gIIGABAB = 0.0 /nIcells
const gICGABAB = 1.5 /nIcells
const count_fMIN = 75
const count_f = 75

frequencies = [1000/0.1,1000/1, 1000/5, 1000/10,1000/15, 1000/20, 1000/25, 1000/30, 1000/35, 1000/40,1000/45, 1000/50]
# frequencies = [ 1000/10,1000/20,1000/40, 1000/50]

# frequencies = [ 50., 33., 25., 20.]
w_init = 0.5
w_Kernel = zeros(length(frequencies))
w_Kernel2 = zeros(length(frequencies))
cd("C:/Users/carol.000/Documents/PhD/Data/SJO Kernels/results_Triplet_SB_H")

for i = 1:length(frequencies)*2

    if i <= length(frequencies)
        freq = frequencies[i]
        delay_pre = 10.

        @time (w_) = simulate_SJO_Triplet(freq, delay_pre, w_init, SB)

        w_Kernel[i] = (w_ / w_init) #*100 / w_[1]
    else
        freq = frequencies[i-length(frequencies)]
        delay_pre = -10.

        @time (w_) = simulate_SJO_Triplet(freq, delay_pre, w_init, SB)

        w_Kernel2[i-length(frequencies)] = (w_ / w_init) #*100 / w_[1]
    end
end

writedlm("wKernel.dat", w_Kernel, header = false)
writedlm("wKernel2.dat", w_Kernel2, header = false)

plot( [0.1,1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], w_Kernel[:] ,color = "blue", legend = false)
plot!([0.1,1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], w_Kernel2[:] ,color = "blue", legend = false)

# plot( [10, 20, 40, 50], w_Kernel[:] ,color = "blue", legend = false)
# plot!([10,20, 40,50], w_Kernel2[:] ,color = "blue", legend = false)
