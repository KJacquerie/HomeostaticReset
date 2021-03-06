
# Julia (.jl) files

All computational experiments are launched thanks to a *Simu_... .jl* calling a function with the differential equations of the conductance based models (neuronal model) and the synaptic weight change (plasticity rules) found in *model_... .jl*. The common methodology/parameters for all codes are explained below:

## Simu_ ... .jl (main code) 

#### Simulation parameters
* *T* is the total duration of the simulation
* *dt* is the time step used for the resolution of ODEs with explicit Euler's method

#### Model parameters (global)
* gather the global parameters used for conductance-based modeling of all cells (Nernst Potentials, membrane capacitances)

#### Model parameters (mean) - Ecells
* gather the parameters used for the ionic conductance of  each cells (E is the presynaptic cell, C the postsynaptic and I the inhibitory cell)

#### Plasticity parameters
* gather the parameters specific of each plasticity rules. The different set of parameters are commented. If the user needs to change the set of parameter of one model, it needs to be done by hand. 
* *wMax* & *wMin* are the hard bounds set for the synaptic weight.

#### Applied currents
* *Iapp* is the constant current applied to the cell, *Istep* a step current applied for a pairing protocol.
* *duration* is the duration of the step current occuring at a given frequency for a pairing protocol.
* *delay_pre* is the $ \Delta t$ imposed between the pre and postsynaptic cell for a pairing protocol.

#### Connectivity
* gather connectivity parameters between the cells such as *AMPA* and *GABA_A*, *GABA_B* connectivities. In this paper, the postsynaptic cell (C) receives *AMPA* input and both pre- and postsynaptic cells have *GABA_A* and *GABA_B* channels activated by the inhibitory cell.

## model_ ... . jl (ODEs)

#### Functions for the differential equations of...
* The membrane voltage, $V_m$, 

$$ C_{m} \frac{V_{m}}{dt} = -\sum I + I_{app} $$
 
where $C_m$ is the membrane capacitance, $I$ is the $i$-th current due to ionic channels $i$, $I_{app}$ is an external applied current. Ionic currents are voltage-dependent; they are expressed as follow:

$$ I = g m^{p}(V_m) h^{q}(V_m) (V_m - E)  $$

where $g$ is the maximal conductance, $m$ is the activation variable, $h$ is the inactivation variable,
$p$ is an integer between 1 and 4, $q$ is either 0 or 1, and $E$ is the reversal potential of the channel.

* The activation and inactivation variables m and p: 

$$  \tau_{m}(V_m)\frac{dm}{dt} = m_{i,\infty}(V_m) - m $$
$$ \tau_{h}(V_m)\frac{dh}{dt} = h_{i,\infty}(V_m) - h $$

* The AMPA and GABA connectivities $ \frac{s_{AMPA}}{dt} $, $ \frac{s_{GABA_A}}{dt} $, $ \frac{s_{GABA_B}}{dt} $. 

* The plasticity rule used

$$ \frac{dw}{dt} = ... $$

(see Methods and Supplementary Materials for more information)

#### function_ ...

A function simulating the computational experiment and simulating the ODEs above with Euler explicit method for each time step *z* and each cell *j*. 
Note: j=1-> presynaptic cell ; j = 2-> inhibitory cell ; j = 3-> postsynaptic cell.

## Data acquisition

* *Simu_scenario* and *model_scenario* used to simulate the scenario of **Figure 1B** with Triplet model.
* *Simu_STDP_..* and *model_STDP_..* used for **Figure 2C** (reproducing the experimental pairing protocol of Bi and Poo, 1998)
* *Simu_SJO_..* and *model_SJO_..* used for **Figure 2C** (reproducing the experimental pairing protocol of Sj??str??m et. al., 2001)
* *Simu_BurstMAT_..* and *model_BurstMAT_..* for **Figure 3A** (simulating the homeostatic reset and saturation for each plasticity rule. Note that for phenomenological models (pair-based and triplet) hard bounds and soft bounds are applied with parameter SB = 0 and SB = 1. )
* *Simu_Variability_..* and *model_BurstMAT_..* for **Figure 3B** (simulating 10 different networks with different ionic conductances (10% variability) for each cell.) 
* *Simu_network* and *model_network* for **Figure 3C** (simulating a 200cells feedforward network with 100 presynaptic cells connected to 100 postsynaptic cells governed by the same inhibitory cell. All initial synaptic weights between pre-post cells are randomly initialized and cells show 20% variability in their $g_{CaT}$ and $g_{KCa} conductances.)

## Models used

Except the scenario codes for **Figure 1B** and network codes for **Figure 3C**, all .jl files are sorted by their model types : 

#### Phenomenological models 
* Pair-based model
* Triplet models 

#### Calcium models 
* **model 1** : Calcium-based synaptic rule with two-triggered thresholds of potentiation and depression. (*Graupner and Brunel, 2016*)
* **model 2 & 3**: *Graupner and Brunel, 2016* model with calcium-dependant time-constant and steady-state values with the set of parameters of the cortex (**model 2**) and the hippocampus (**model 3**).
* **model 4 & 5**: *Graupner and Brunel, 2012* model with the set of parameters of the cortex (**model 4**) and the hippocampus (**model 5**). 
* **model 6**: Calcium-dependent plasticity rule with a continuous function of the calcium concentration for the cortex. (*Shouval, Bear and Cooper, 2002*). 
* **model 7**: Calcium-dependent plasticity rule with coupling between the synaptic weight and the calcium influx. (*Deperrois and Graupner, 2020*)
* **model 8**:  Calcium-dependent plasticity rule with coupling between the synaptic weight and the calcium influx and a short-term depression at the presynaptic site (*Deperrois and Graupner, 2020*).
* **model 9, 10, 11, 12, 13, 14**: Calcium-based plasticity rule with hard bounds. *Graupner and Brunel, 2016* (**model 9**), *Shouval, Bear and Cooper, 2002* (**model 10**), *Deperrois and Graupner, 2020* with coupling between the synaptic weight and the calcium influx (Threshold **model 11**, continuous form **model 12**),  *Deperrois and Graupner, 2020* with coupling between the synaptic weight and the calcium influx and a short-term depression at the presynaptic site  (Threshold **model 13**, continuous form **model 14**)
