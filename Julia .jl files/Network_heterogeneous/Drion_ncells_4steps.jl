# gating functions
boltz(V::Float64,A::Float64,B::Float64) = 1/(1 + exp((V+A)/B))
tauX(V::Float64,A::Float64,B::Float64,D::Float64,E::Float64) = A - B/(1+exp((V+D)/E))
mNainf(V::Float64) = boltz(V,35.5,-5.29)
taumNa(V::Float64) = tauX(V,1.32,1.26,120.,-25.)
hNainf(V::Float64) = boltz(V,48.9,5.18)
tauhNa(V::Float64) = (0.67/(1+exp((V+62.9)/-10.0)))*(1.5 + 1/(1+exp((V+34.9)/3.6)))
mKdinf(V::Float64) = boltz(V,12.3,-11.8)
taumKd(V::Float64) = tauX(V,7.2,6.4,28.3,-19.2)
mCaTinf(V::Float64) = boltz(V,67.1,-7.2)
taumCaT(V::Float64) = tauX(V,21.7,21.3,68.1,-20.5)
hCaTinf(V::Float64) = boltz(V,80.1,5.5)
tauhCaT(V::Float64) = 2*tauX(V,205.,89.8,55.,-16.9)
mHinf(V::Float64) = boltz(V,80.,6.)
taumH(V::Float64) = tauX(V,272.,-1149.,42.2,-8.73)
mKCainf(Ca::Float64) = (Ca/(Ca+Kd))^2
Tm(V::Float64) = 1/(1+exp(-(V-2)/5))

ICaT(V::Float64,mCaT::Float64,hCaT::Float64,gCaT::Float64) = gCaT*mCaT^3*hCaT*(V-VCa)

function dV(gNa::Float64, gKd::Float64, gl::Float64, V::Float64, mNa::Float64, hNa::Float64, mKd::Float64, mCaT::Float64, hCaT::Float64, mH::Float64, Ca::Float64, Iapp::Float64, Istep::Float64, gCaT::Float64, gH::Float64, gKCa::Float64)
  (dt)*(1/C)*(-gNa*mNa^3*hNa*(V-VNa) -gKd*mKd^4*(V-VK) -gl*(V-Vl) -ICaT(V,mCaT,hCaT,gCaT) -gKCa*mKCainf(Ca)*(V-VK) -gH*mH*(V-VH) +Iapp +Istep)
end

dmNa(V::Float64,mNa::Float64) = (dt)*((1/taumNa(V))*(mNainf(V) - mNa))
dhNa(V::Float64,hNa::Float64) = (dt)*((1/tauhNa(V))*(hNainf(V) - hNa))
dmKd(V::Float64,mKd::Float64) = (dt)*((1/taumKd(V))*(mKdinf(V) - mKd))
dmCaT(V::Float64,mCaT::Float64) = (dt)*((1/taumCaT(V))*(mCaTinf(V) - mCaT))
dhCaT(V::Float64,hCaT::Float64) = (dt)*((1/tauhCaT(V))*(hCaTinf(V) - hCaT))
dmH(V::Float64,mH::Float64) = (dt)*((1/taumH(V))*(mHinf(V) - mH))
dCa(V::Float64,mCaT::Float64,hCaT::Float64,Ca::Float64,gCaT::Float64,k1::Float64,k2::Float64) = (dt)*(-k1*ICaT(V,mCaT,hCaT,gCaT) -k2*Ca)

dAMPA(V::Float64,AMPA::Float64) = (dt)*(1.1*Tm(V)*(1-AMPA)-0.19*AMPA)
dGABAA(V::Float64,GABAA::Float64) = (dt)*(0.53*Tm(V)*(1-GABAA)-0.18*GABAA)
dGABAB(V::Float64,GABAB::Float64) = (dt)*(0.016*Tm(V)*(1-GABAB)-0.0047*GABAB)


#=
function dx(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_p + 1*(Vprev<= 0 && V >0)
end

function dx2(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_x + 1*(Vprev<= 0 && V >0)
end

function dy(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_m + 1*(Vprev<= 0 && V >0)
end

function dy2(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_y + 1*(Vprev<= 0 && V >0)
end


function dw(new::Float64, w::Float64, tau::Float64)
  (dt)*(new - w)*1/tau
end
=#

function simulateTOY_ncells(ncells::Int64,nPrecells::Int64,nPostcells::Int64,IappE::Float64,IappI::Float64,TstepEinit::Int64,TstepEfinal::Int64,IstepE::Float64,TstepIinit1::Int64,TstepIinit2::Int64,TstepIinit3::Int64,TstepIinit4::Int64,TstepIfinal::Int64,IstepI1::Float64,IstepI2::Float64,IstepI3::Float64,IstepI4::Float64,gAMPA::Float64, gGABAA, gGABAB, wmat)

  # Initial conditions
  println(wmat)

  V=-60*ones(ncells)
  Vprev=-60*ones(ncells)
  mNa=mNainf(V[1])*ones(ncells)
  hNa=hNainf(V[1])*ones(ncells)
  mKd=mKdinf(V[1])*ones(ncells)
  mCaT=mCaTinf(V[1])*ones(ncells)
  hCaT=hCaTinf(V[1])*ones(ncells)
  mH=mHinf(V[1])*ones(ncells)
  Ca=((-k1_E/k2_E)*ICaT(V[1],mCaT[1],hCaT[1],gCaT_E[1]))*ones(ncells)
  AMPA=zeros(ncells)
  GABAA=zeros(ncells)
  GABAB=zeros(ncells)

##Plasticity parameters
  #=
  x = zeros(nPrecells)
  x_prev = zeros(nPrecells)
  y = zeros(nPostcells)
  y_prev = zeros(nPostcells)
  x2 = zeros(nPrecells)
  x2_prev = zeros(nPrecells)
  y2 = zeros(nPostcells)
  y2_prev = zeros(nPostcells)
  pre_spikes = zeros(nPrecells)
  post_spikes = zeros(nPostcells)
  =#

  cpre_var =zeros(nPrecells)
  cpre_prev=zeros(nPrecells)
  cpost_var = zeros(nPostcells)
  cpost_prev = zeros(nPostcells)

  event_pre = zeros(nPrecells)
  event_post = zeros(nPostcells)

  ca_var =zeros(nPostcells, nPrecells)

##


  VV = zeros(Tdt,ncells)
  ww = zeros(Tdt,length(wmat))

  Vold = zeros(convert(Int64, floor(D_pre/dt+1)),nPrecells)

  # Step start and stop values
  TstartE::Int64 = convert(Int64,TstepEinit/dt)
  TstopE::Int64 = convert(Int64,TstepEfinal/dt)
  TstartC::Int64 = convert(Int64,tstepCinit/dt)
  TstopC::Int64 = convert(Int64,tstepCfinal/dt)

  TstartI1::Int64 = convert(Int64,TstepIinit1/dt)
  TstartI2::Int64 = convert(Int64,TstepIinit2/dt)
  TstartI3::Int64 = convert(Int64,TstepIinit3/dt)
  TstartI4::Int64 = convert(Int64,TstepIinit4/dt)
  TstopI::Int64   = convert(Int64,TstepIfinal/dt)

  for z = 1:Tdt
    #Isyn=zeros(ncells)

    for j = 1:ncells

##Pre cells (E)
      if j<=nPrecells
        Iapp = IappE
        if z >= TstartE && z<= TstopE
          Iappstep = IstepE
        else
          Iappstep = 0.
        end
        V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])
      end


## Post cells (C)
      if j>nPrecells && j<ncells
        Iapp = IappC
        if z >= TstartC && z<= TstopC
          Iappstep = IstepC
        else
          Iappstep = 0.
        end
        V[j] += dV(gNavec_C[j], gKdvec_C[j], glvec_C[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_C[j], gHvec_C[j], gKCavec_C[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_C[j],k1vec_C[j], k2vec_C[j])
      end


## I Cells (I)

      if j==ncells
        Iapp = IappI
        if z >= TstartI1 && z< TstartI2
          Iappstep = IstepI1
        elseif z >= TstartI2 && z< TstartI3
          Iappstep = IstepI2
        elseif z >= TstartI3 && z< TstartI4
          Iappstep = IstepI3
        elseif z >= TstartI4 && z< TstopI
          Iappstep = IstepI4
        else
          Iappstep = 0.
        end
        V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
      end

## Connectivities
      if j <=nPrecells
            ##
            #Isyn[j] += gGABAA[j]*GABAA[ncells]*(Vprev[j]+70)+gGABAB[j]*GABAB[ncells]*(Vprev[j]+85)
            V[j] += (dt)*(1/C)*(-gGABAA[j]*GABAA[ncells]*(Vprev[j]+70))
            V[j] += (dt)*(1/C)*(-gGABAB[j]*GABAB[ncells]*(Vprev[j]+85))

            if(z> (D_pre/dt+1) )
              if(Vold[end-1,j]>=0 && Vold[end,j]<0)
                event_pre[j]=1.
              else
                event_pre[j]=0.
              end
            end
            Vold[2:end,j] = Vold[1:end-1,j]
            Vold[1,j] = copy(V[j])

            cpre_var[j] +=(dt)*(1)*(-cpre_prev[j]/tau_Ca) +event_pre[j]*C_Pre


      end

      if j > nPrecells && j<ncells

          #Isyn[j] += gGABAA[j]*GABAA[ncells]*(Vprev[j]+70)+gGABAB[j]*GABAB[ncells]*(Vprev[j]+85)
          V[j] += (dt)*(1/C)*(-gGABAA[j]*GABAA[ncells]*(Vprev[j]+70))
          V[j] += (dt)*(1/C)*(-gGABAB[j]*GABAB[ncells]*(Vprev[j]+85))
          for k=1:nPrecells
              #Isyn[j]+= gAMPA*wmat[j-nPrecells,k]*AMPA[k]*(Vprev[j]-0)
              V[j] += (dt)*(1/C)*(-gAMPA*wmat[j-nPrecells,k]*AMPA[k]*(Vprev[j]-0))
          end

          if(V[j]>=0. && Vprev[j]<0.)
            event_post[j-nPrecells]=1.
            #println(event_post)
          else
            event_post[j-nPrecells]=0.
          end
          cpost_var[j-nPrecells] += (dt)*(1)*(-cpost_prev[j-nPrecells]/tau_Ca)+event_post[j-nPrecells]*C_Post


      end

## Plasticity
  #=
      if (j <=nPrecells)

        pre_spikes[j]=0
        x_prev[j] = x[j]
        x[j] += dx(x[j], V[j], Vprev[j])
        x2_prev[j] = x2[j]
        x2[j] += dx2(x2[j], V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          pre_spikes[j] = 1  #Pré en dernier
        end

      end
      =#

      #=
      if(j>nPrecells && j< ncells)

        post_spikes[j-nPrecells]=0

        y_prev[j-nPrecells] = y[j-nPrecells]
        y[j-nPrecells] += dy(y[j-nPrecells], V[j], Vprev[j])
        y2_prev[j-nPrecells] = y2[j-nPrecells]
        y2[j-nPrecells] += dy2(y2[j-nPrecells], V[j], Vprev[j])
        if (V[j] >= 0 && Vprev[j]< 0)
          post_spikes[j-nPrecells] = 1 #Post en dernier
        end
      end
      =#


      if (j>nPrecells && j< ncells)

        if(z>TstartI1)
          #println("here")
          for k=1:nPrecells
            if(ca_var[j-nPrecells, k]>= theta_p)
              wmat[j-nPrecells, k] +=(dt)*(1/tauw_p)*(Omega_p - wmat[j-nPrecells, k])
            elseif(ca_var[j-nPrecells, k]>=theta_d && ca_var[j-nPrecells, k]<theta_p)
              wmat[j-nPrecells, k] +=(dt)*(1/tauw_d)*(Omega_d - wmat[j-nPrecells, k])
            else
              wmat[j-nPrecells, k] +=(dt)*0.
            end
            ca_var[j-nPrecells, k] = cpre_prev[k] + cpost_prev[j-nPrecells]
            cpre_prev[k] = cpre_var[k]
          end
        end
        cpost_prev[j-nPrecells] = cpost_var[j-nPrecells]
      end


      #=
      if (j>nPrecells && j< ncells)
        for k=1:nPrecells
          if pre_spikes[k] == 1
            if SB == 0
              new = -y_prev[j-nPrecells]*A2_m  + wmat[j-nPrecells,k]
            else
              new = -y_prev[j-nPrecells]*wmat[j-nPrecells,k]*A2_m  + wmat[j-nPrecells,k]
            end
            wmat[j-nPrecells,k] = new
            if (wmat[j-nPrecells,k] > wMax)
              wmat[j-nPrecells,k]= wMax
            end
            if (wmat[j-nPrecells,k] < wMin)
              wmat[j-nPrecells,k]= wMin
            end

          end
        end
      end
      =#
      #=
      if (j<=nPrecells)
        for k=1:nPostcells
          if post_spikes[k] == 1
            if SB == 0
              new = x_prev[j]*(A2_p+ A3_p*y2_prev[j]) + wmat[k,j]
            else
              new = x_prev[j]*(1-wmat[k,j])*(A2_p+ A3_p*y2_prev[j]) + wmat[k,j]
            end
            wmat[k,j] = new
            if (wmat[k,j] > wMax)
              wmat[k,j]= wMax
            end
            if (wmat[k,j] < wMin)
              wmat[k,j]= wMin
            end

          end
        end
      end
      =#


##



      mNa[j] += dmNa(Vprev[j],mNa[j])
      hNa[j] += dhNa(Vprev[j],hNa[j])
      mKd[j] += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j] += dmH(Vprev[j],mH[j])

      AMPA[j] += dAMPA(Vprev[j],AMPA[j])
      GABAA[j] += dGABAA(Vprev[j],GABAA[j])
      GABAB[j] += dGABAB(Vprev[j],GABAB[j])

      Vprev = copy(V)


    end

    VV[z,:] = copy(V')

    for i=1:length(wmat)
      ww[z,i] = wmat[i]
    end

  end
  return VV, ww
end
