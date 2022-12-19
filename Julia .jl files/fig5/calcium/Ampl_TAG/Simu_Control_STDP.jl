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
include("model_Control_STDP.jl")

# Simulation parameters

const T = 1000000
const dt = 0.01
const Tdt = convert(Int64, T / dt)
const Ttransient = convert(Int64, 1 / dt)
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



## Param


const expm = "cortex"
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

#=
const expm="hippocampus"
const gEC_const = 0.01
const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410
const C_Post = 1.62138
const D_pre  = 9.53709 #[ms]

const tau_w    = 520.76129e3 #[s>ms]
const gamma_p  = 597.08922
const gamma_d  = 137.7586
const n_lin = 1.
const theta_p = 1.45   #2.009289
const theta_d = 1.
=#


#Omega_p = gamma_p/(gamma_p+gamma_d)
Omega_d = 0.
Omega_0 = 0.

tauw_p = tau_w /(gamma_p + gamma_d)
tauw_d = tau_w / gamma_d
tauw_0 = 0.




##
const count_f = 60
freq_intensity = 1
freq = 1000/freq_intensity

delays_mat = range(-50., 50., length = 51)



w_init = [ 0.1 0.5 0.9]
w_name =[ 1 5 9]

wfinal = zeros(length(delays_mat), length(w_init))
wdelta_ = zeros(length(delays_mat), length(w_init))

for idx_w = 1:length(w_init)
    # RUN MODEL
    for idx_loop=1:length(delays_mat)
        delay_pre=delays_mat[idx_loop]
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
            delay_pre,
            w_init[idx_w]
        )
        wfinal[idx_loop,idx_w] = w[idx_stop]
    end

    wdelta_[:,idx_w] = wfinal[:,idx_w]/w_init[idx_w]
end

writedlm(@sprintf("/Users/kathleen/Documents/PhD/2022-Paris/Julia/omega_winit/results_STDP/%s/STDP_omega.dat", expm), wdelta_, header=false)

plot(delays_mat,wdelta_[:,1]*w_init[1] - w_init[1]*ones(size(wdelta_,1)) , color = "blue", xlabel = "delta t (ms)", ylabel = " delta w (-)", legend= (frameon=false))
plot!(delays_mat,wdelta_[:,2]*w_init[2] - w_init[2]*ones(size(wdelta_,1)), color = "black", xlabel = "delta t (ms)", ylabel = " delta w (-)", legend= (frameon=false))
plot!(delays_mat,wdelta_[:,3]*w_init[3] - w_init[3]*ones(size(wdelta_,1)), color = "red", xlabel = "delta t (ms)", ylabel = " delta w (-)", legend= (frameon=false))
