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
include("model_CC_rel_BurstMATD.jl")
include("Drion_PARAMS_4steps.jl")

# Simulation parameters

const T = 80000
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
const IstepE2 = 0.
const tstepEinit = 500 # Transient
const tstepEfinal = T

##Cellule I
const IappI = 3.
#const IstepI1 = -0.9-IappI
const IstepI2 = 0.
const IstepI3 = 0.
const IstepI4 = 0.
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

const tau_Ca = 22.27212 #[ms]
const C_Pre  = 0.84410
const C_Post = 1.62138
const D_pre  = 9.53709 #[ms]

#const tau_w    = 520.76129e3 #[s>ms]
const gamma_p  = 597.08922
const gamma_d  = 137.7586
const n_lin = 1.
const theta_p = 2.009289
const theta_d = 1.

const eta_lin = ((n_lin*(C_Pre+C_Post)-C_Post)/C_Pre )-1


## Parameters for Shouval Rule

const expm="CC_rel"
const strength = "weak"
if(strength=="weak")
    const gEC_const = 0.01
else
    const gEC_const = 0.4
end
a1 = 1.
a2 = 2.
m1 = 0.25
b1 = 40
b2 = 10#40
m2=0.5

P1=0.1*4e4
P2=P1*1e-6
P3=2.4
P4=1

wMAX = 1.



##

const freq_intensity=10.
const freq = 1000/freq_intensity
const delay_pre = 10.0#-D_pre#-10.0-D_pre


# RUN MODEL
const n_net = 1 # one circuit atm
const j_net = 1 #just to tell the idx of the ntk we are considering
const w_init_vec=0.0:0.1:1.0
const IstepI_vec = -1.7-IappI:0.1:-0.9-IappI

function simu_wI()
    for j=1:1:length(w_init_vec)
        w_init=w_init_vec[j]
        println("wi=",w_init)
        for p=1:1:length(IstepI_vec)
            IstepI1 = IstepI_vec[p]
            cd(@sprintf("/Users/kathleen/Documents/PhD/2022-reset/model10/results_sleep/%s/w%d/%d/",strength,w_init_vec[j]*10,p))


            #Initialisation des vecteurs
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


            @time (Vconnect_spk, w, wf) = simulateTOY_ncells(
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
                w_init
            )

            (ISIs_depol, ISIs_hyperpol, PARAMS_depol, PARAMS_hyperpol, freq_depol, freq_hyperpol) = compute_params(Vconnect_spk[:,:],Ttransient,tstepIinit1,tstepIinit1,tstepIinit2)
            SPB_depol1[:,j_net] = PARAMS_depol[:,1]
            SPB_hyperpol1[:,j_net] = PARAMS_hyperpol[:,1]
            PER_depol1[:,j_net] = PARAMS_depol[:,2]
            PER_hyperpol1[:,j_net] = PARAMS_hyperpol[:,2]
            DC_depol1[:,j_net] = PARAMS_depol[:,3]
            DC_hyperpol1[:,j_net] = PARAMS_hyperpol[:,3]
            IBF_depol1[:,j_net] = PARAMS_depol[:,4]
            IBF_hyperpol1[:,j_net] = PARAMS_hyperpol[:,4]
            freq_depol_vec1[:,j_net] = freq_depol
            freq_hyperpol_vec1[:,j_net] = freq_hyperpol

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

           writedlm("Vspk.dat", Vconnect_spk, header = false)
           writedlm("w.dat", w, header = false)
           writedlm("wf.dat", wf, header = false)


        end # loop on IstepI

    end # loop on w_init
end


@time  simu_wI()
