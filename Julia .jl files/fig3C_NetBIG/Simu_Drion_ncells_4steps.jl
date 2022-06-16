# Defines output directory

# Loads packages
using Plots
using DelimitedFiles
using Statistics
using Images, ImageView

# Include model
include("Drion_ncells_4steps.jl")
include("Drion_PARAMS_4steps.jl")


# Simulation parameters
const T = 50000
const dt = 0.01
const Tdt = convert(Int64,T/dt)
const Ttransient= convert(Int64,400/dt)
const t = range(dt,T,length=Tdt)

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
const gNa_E = 170.
const gKd_E = 40.
const k1_E = 1.e-1
const k2_E = 0.1e-1
const gH_E = 0.01
const gKCa_E = 4.
const gCaT_E = 0.55

# Model parameters (mean) - Icells
const gl_I = 0.055
const gNa_I = 170.
const gKd_I = 40.
const k1_I = 1.e-1
const k2_I = 0.1e-1
const gH_I = 0.01
const gKCa_I = 4.
const gCaT_I = 0.55

# Model parameters (mean) - Ecells
const gl_C = 0.055
const gNa_C = 170.
const gKd_C = 40.
const k1_C = 1.e-1
const k2_C = 0.1e-1
const gH_C = 0.01
const gKCa_C = 4.
const gCaT_C = 0.55


# Simulations
const nPrecells = 3  # Number of excitatory cells
const nIcells = 1    # Number of inhibitory cells
const nPostcells = 3 # Number of excitatory cells
const ncells = nPrecells+nIcells+nPostcells

variab = 0.2

const k1vec_E = k1_E*ones(ncells)
const k2vec_E = k2_E*ones(ncells)

const k1vec_I = k1_I*ones(ncells)
const k2vec_I = k2_I*ones(ncells)

const k1vec_C = k1_C*ones(ncells)
const k2vec_C = k2_C*ones(ncells)

gNavec_E = gNa_E*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 70*(rand(ncells)-0.5)
gKdvec_E = gKd_E *ones(ncells)#* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 40*(rand(ncells)-0.5)
gHvec_E = gH_E*ones(ncells)# * (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.002*(rand(ncells)-0.5)
gCaTvec_E = gCaT_E .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))  #(1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.35*(rand(ncells)-0.5)
gKCavec_E = gKCa_E .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))  #(1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 2*(rand(ncells)-0.5)
glvec_E = gl_E*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.015*(rand(ncells)-0.5)

gNavec_I = gNa_I*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 70*(rand(ncells)-0.5)
gKdvec_I = gKd_I*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 40*(rand(ncells)-0.5)
gHvec_I = gH_I*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 0.002*(rand(ncells)-0.5)
gCaTvec_I = gCaT_I.*ones(ncells) #.* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 0.35*(rand(ncells)-0.5)
gKCavec_I = gKCa_I.*ones(ncells) #.* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 2*(rand(ncells)-0.5)
glvec_I = gl_I*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.015*(rand(ncells)-0.5)

gNavec_C = gNa_C*ones(ncells) #*(1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 70*(rand(ncells)-0.5)
gKdvec_C = gKd_C*ones(ncells)# *(1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 40*(rand(ncells)-0.5)
gHvec_C = gH_C*ones(ncells)# * (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.002*(rand(ncells)-0.5)
gCaTvec_C = gCaT_C .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5)) #- 0.35*(rand(ncells)-0.5)
gKCavec_C = gKCa_C .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 2*(rand(ncells)-0.5)
glvec_C = gl_C*ones(ncells) #* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))#- 0.015*(rand(ncells)-0.5)

#simulateTOY_ncells(ncells::Int64,nEcells::Int64,nIcells::Int64,Iapp::Float64,Tstepinit::Int64,Tstepfinal::Int64,Istep::Float64,gEE::Float64,gEI::Float64,gIE::Float64,gII::Float64,gIE2::Float64,gII2::Float64)
const IappE = 0.
const IstepE = 0.
const tstepEinit = T
const tstepEfinal = T

const IappC = 0.
const IstepC = 0.
const tstepCinit = T
const tstepCfinal = T

const IappI = 3.
const IstepI1 = - IappI - 1.2
const IstepI2 = 0.
const IstepI3 = 0.
const IstepI4 = 0.
const tstepIinit1 = 400
const tstepIinit2 = T
const tstepIinit3 = T
const tstepIinit4 = T
const tstepIfinal = copy(T)

## Plasticity parameters
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

#= #Minimal model Pfister2006 Hippocampe (Hard bounds)
const A2_m = 3.5e-3#7e-3
const A3_m = 0.#2.3e-4 #0
const A2_p = 5.3e-3#5e-10 #0
const A3_p = 8e-3#6.5e-3

const tau_m = 33.7
const tau_p = 16.8
const tau_x = 101#101#101
const tau_y = 40.#125
const SB = 1
=#

#=
const wMax = 1.
const wMin = 0.
=#
## connectivity

const gAMPA = 0.01/nPrecells

gGABAA_vec = 1.5/nIcells .* [0.8, 1.1, 1.4, 1, 1, 1, 1] #ones(ncells)# .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))
gGABAB_vec = 2.0/nIcells .* [0.8, 1.1, 1.4, 1, 1, 1, 1]# .* (1 .-2*variab*(rand(Float64, (1, ncells)).-0.5))

w_mat = rand(nPostcells, nPrecells)

# RUN MODEL
println([gCaTvec_E[1:nPrecells]' gCaTvec_C[nPrecells+1:ncells-1]' gCaTvec_I[ncells];gKCavec_E[1:nPrecells]' gKCavec_C[nPrecells+1:ncells-1]' gKCavec_I[ncells]; gGABAA_vec'; gGABAB_vec'])


#Initialisation des vecteurs
n_net = 1
SPB_depol1 = zeros(ncells,n_net)
SPB_hyperpol1 = zeros(ncells,n_net)
PER_depol1 = zeros(ncells,n_net)
PER_hyperpol1 = zeros(ncells,n_net)
DC_depol1 = zeros(ncells,n_net)
DC_hyperpol1 = zeros(ncells,n_net)
IBF_depol1 = zeros(ncells,n_net)
IBF_hyperpol1 = zeros(ncells,n_net)
freq_depol_vec1 = zeros(ncells,n_net)
freq_hyperpol_vec1 = zeros(ncells,n_net)

@time (Vconnect_spk, w_) = simulateTOY_ncells(ncells,nPrecells,nPostcells,IappE,IappI,tstepEinit,tstepEfinal,IstepE,tstepIinit1,tstepIinit2,tstepIinit3,tstepIinit4,tstepIfinal,IstepI1,IstepI2,IstepI3,IstepI4,gAMPA, gGABAA_vec, gGABAB_vec, w_mat)

nEcells=nPrecells
(ISIs_depol, ISIs_hyperpol, PARAMS_depol, PARAMS_hyperpol, freq_depol, freq_hyperpol) = compute_params(Vconnect_spk,Ttransient,tstepIinit1,tstepIinit1,tstepIinit2)
SPB_depol1[:,1]     =   PARAMS_depol[:,1]
SPB_hyperpol1[:,1]  = PARAMS_hyperpol[:,1]
PER_depol1[:,1]     = PARAMS_depol[:,2]
PER_hyperpol1[:,1]  = PARAMS_hyperpol[:,2]
DC_depol1[:,1]      = PARAMS_depol[:,3]
DC_hyperpol1[:,1]   = PARAMS_hyperpol[:,3]
IBF_depol1[:,1]     = PARAMS_depol[:,4]
IBF_hyperpol1[:,1]  = PARAMS_hyperpol[:,4]
freq_depol_vec1[:,1]    = freq_depol
freq_hyperpol_vec1[:,1] = freq_hyperpol

writedlm("SPB_depol1.dat", SPB_depol1, header = false)
writedlm("SPB_hyperpol1.dat", SPB_hyperpol1, header = false)
writedlm("PER_depol1.dat", PER_depol1, header = false)
writedlm("PER_hyperpol1.dat", PER_hyperpol1, header = false)
writedlm("DC_depol1.dat", DC_depol1, header = false)
writedlm("DC_hyperpol1.dat", DC_hyperpol1, header = false)
writedlm("IBF_depol1.dat", IBF_depol1, header = false)
writedlm("IBF_hyperpol1.dat", IBF_hyperpol1, header = false)
writedlm("freq_depol1.dat", freq_depol_vec1, header = false)
writedlm("freq_hyperpol1.dat", freq_hyperpol_vec1, header = false)



## Uncomment to save data
# cd("C:/Users/...")
# writedlm("w.dat", w_, header = false)
# writedlm("wf.dat", w_[end,:], header = false)
# writedlm("winit.dat", w_[1,:], header = false)
# writedlm("gVar.dat", [gCaTvec_E[1:nPrecells]' gCaTvec_C[nPrecells+1:ncells-1]' gCaTvec_I[ncells];gKCavec_E[1:nPrecells]' gKCavec_C[nPrecells+1:ncells-1]' gKCavec_I[ncells]; gGABAA_vec; gGABAB_vec], header = false)


## Plot V
#=
Vplot1 = plot(t[100000:350000],Vconnect[100000:350000,1], color="blue")
Vplot2 = plot(t[100000:350000],Vconnect[100000:350000,2], color="blue")
Vplot3 = plot(t[100000:350000],Vconnect[100000:350000,3], color="blue")

Vplot4 = plot(t[100000:350000],Vconnect[100000:350000,4], color="grey")
Vplot5 = plot(t[100000:350000],Vconnect[100000:350000,5], color="grey")
Vplot6 = plot(t[100000:350000],Vconnect[100000:350000,6], color="grey")

Vplot7 = plot(t[100000:350000],Vconnect[100000:350000,7], color="red")

plot(Vplot1,Vplot4, Vplot7,Vplot2, Vplot5, Vplot7, Vplot3,  Vplot6,  Vplot7, layout=(3,3))
=#
## Plot w#

plot(t, w_[:,1])

##

writedlm("/Users/kathleen/Documents/PhD/2022-NMOD/Network_heterogeneousBIG/trace/Vconnec_spk.dat", Vconnect, header = false)
writedlm("/Users/kathleen/Documents/PhD/2022-NMOD/Network_heterogeneousBIG/trace/w.dat", w_, header = false)
writedlm("/Users/kathleen/Documents/PhD/2022-NMOD/Network_heterogeneousBIG/trace/gCaTvec_E.dat", gCaTvec_E, header = false)
writedlm("/Users/kathleen/Documents/PhD/2022-NMOD/Network_heterogeneousBIG/trace/gCaTvec_C.dat", gCaTvec_C, header = false)
writedlm("/Users/kathleen/Documents/PhD/2022-NMOD/Network_heterogeneousBIG/trace/gCaTvec_I.dat", gCaTvec_I, header = false)
