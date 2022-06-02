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

function Heaviside_p(ca::Float64, threshold_p::Float64)
  if(ca >=threshold_p)
    res=1.
  else
    res=0.
  end
  return res
end

function Heaviside_d(ca::Float64, threshold_d::Float64)
  if ca >= threshold_d #&& ca < threshold_p
    res = 1.
  else
    res = 0.
  end
  return res
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

dcpre(cpre::Float64, event_pre::Float64) = (dt)*(-(cpre/tau_Ca)) + event_pre*C_Pre
dcpost(cpost::Float64,cpre::Float64, event_post::Float64) = (dt)*(-(cpost/tau_Ca)) + event_post*C_Post + eta_lin*event_post*cpre


function simulateTOY_ncells(ncells::Int64,nEcells::Int64,nIcells::Int64,IappE::Float64,IappI::Float64,TstepEinit::Int64,TstepEfinal::Int64,IstepE::Float64,TstepIinit1::Int64,TstepIinit2::Int64,TstepIinit3::Int64,TstepIinit4::Int64,TstepIfinal::Int64,IstepI1::Float64,IstepI2::Float64,IstepI3::Float64,IstepI4::Float64, gEC::Float64, freq::Float64, delay_pre::Float64)
  # Initial conditions
  V=-60*ones(ncells)
  Vprev=-60*ones(ncells)
  x::Float64 = 0.
  y::Float64 = 0.
  mNa=mNainf(V[1])*ones(ncells)
  hNa=hNainf(V[1])*ones(ncells)
  mKd=mKdinf(V[1])*ones(ncells)
  mCaT=mCaTinf(V[1])*ones(ncells)
  hCaT=hCaTinf(V[1])*ones(ncells)
  mH=mHinf(V[1])*ones(ncells)
  Ca=((-k1_E/k2_E)*ICaT(V[1],mCaT[1],hCaT[1],gCaT_E[1]))*ones(ncells)
  AMPA=zeros(ncells)
  GABAA = zeros(ncells)
  GABAB = zeros(ncells)
  sSYN = zeros(ncells)

  cpre = zeros(Tdt)
  cpre_var =0.
  cpre_prev=0.
  cpost = zeros(Tdt)
  cpost_var = 0.
  cpost_prev=0.
  ca_var =0.
  c_plot = zeros(Tdt)

  #flag_Vpre=1
  #flag_Vpost=1

  idx_pre = 1
  idx_post= 1

  event_pre=0.
  event_post=0.

  Vold = zeros(convert(Int64, floor(D_pre/dt+1)),1)


  ww=zeros(Tdt)
  ww_var = 0.5
  count_post=0.
  count_pre =0.
  idx_stop = 0.
  #VV = zeros(Tdt,ncells)

  # Step start and stop values
  TstartE::Int64 = convert(Int64,TstepEinit/dt)
  TstopE::Int64 = convert(Int64,TstepEfinal/dt)
  TstartI1::Int64 = convert(Int64,TstepIinit1/dt)
  TstartI2::Int64 = convert(Int64,TstepIinit2/dt)
  TstartI3::Int64 = convert(Int64,TstepIinit3/dt)
  TstartI4::Int64 = convert(Int64,TstepIinit4/dt)
  TstopI::Int64   = convert(Int64,TstepIfinal/dt)

  for z = 1:Tdt

    for j = 1:ncells

      # Excitatory cell
      if (j<=nEcells)
        Iapp = IappE

        z_delayed = z*dt + delay_pre
        if (mod(z_delayed, freq) > freq - duration && z_delayed >=0 )
          Iappstep = IstepE
        else
          Iappstep = 0.
        end

        if z >= TstartE && z< TstopE
          Iappstep = IstepE2
        end

        V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])
      end


      # Inhibitory cell
      if (j > nEcells && j <= nEcells+nIcells)
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

      #  Cortical cell
      if ( j > nIcells+nEcells && j <=ncells)
        Iapp = IappC
        z_delayed = z*dt

        if (mod(z_delayed, freq) > freq - duration && z_delayed >=0 )
          Iappstep = IstepC
        else
          Iappstep = 0.
        end


        if (z >= tStartC/dt && z< tStopC/dt)
          Iappstep = IstepC2
        end

        V[j] += dV(gNavec_I[j], gKdvec_I[j], glvec_I[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_I[j], gHvec_I[j], gKCavec_I[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_I[j],k1vec_I[j], k2vec_I[j])
      end

      # Inhibition
        if j == 1
          V[j] += (dt)*(1/C)*(-gIEGABAA*GABAA[2]*(Vprev[j]+70))
          V[j] += (dt)*(1/C)*(-gIEGABAB*GABAB[2]*(Vprev[j]+85))

          if(V[j]>=0 && Vprev[j]<0)
            count_pre +=1
            #event_pre=1.
          #else
            #event_pre=0.
          end
          if(z> (D_pre/dt+1) )
            if(Vold[end-1]>=0 && Vold[end]<0)
              event_pre=1.
            else
              event_pre=0.
            end
          end
          Vold[2:end] = Vold[1:end-1]
          Vold[1] = copy(V[1])

        end

        if j == 3
          V[j] += (dt)*(1/C)*(-gICGABAA*GABAA[2]*(Vprev[j]+70))
          V[j] += (dt)*(1/C)*(-gICGABAB*GABAB[2]*(Vprev[j]+85))

        end

        # Excitation
        if j == 3
          V[j] += (dt)*(1/C)*(-gEC*ww_var*AMPA[1]*(Vprev[j]-0))

          if(V[j]>=0 && Vprev[j]<0)
            count_post +=1
            event_post=1.
          else
            event_post=0.
          end



        end

      mNa[j]  += dmNa(Vprev[j],mNa[j])
      hNa[j]  += dhNa(Vprev[j],hNa[j])
      mKd[j]  += dmKd(Vprev[j],mKd[j])
      mCaT[j] += dmCaT(Vprev[j],mCaT[j])
      hCaT[j] += dhCaT(Vprev[j],hCaT[j])
      mH[j]   += dmH(Vprev[j],mH[j])

      AMPA[j]  += dAMPA(Vprev[j],AMPA[j])
      GABAA[j] += dGABAA(Vprev[j],GABAA[j])
      GABAB[j] += dGABAB(Vprev[j],GABAB[j])

      Vprev = copy(V)


    end


    cpre_var  += (dt)*(1)*(-cpre_prev/tau_Ca) +event_pre*C_Pre
    cpost_var += (dt)*(1)*(-cpost_prev/tau_Ca)+event_post*C_Post #+ eta_lin*event_post*cpre_prev


    ww_var +=(dt)*(1/tau_w)*(-ww_var*(1-ww_var)*(wfix-ww_var)+gamma_p*(1-ww_var)*Heaviside_p(ca_var,theta_p)-gamma_d*ww_var*Heaviside_d(ca_var, theta_d))


    ca_var = cpre_var +cpost_var
    cpre_prev = cpre_var
    cpost_prev= cpost_var

    ww[z] = copy(ww_var)

    if freq == 10000.
      if count_pre>=count_fMIN &&count_post>=count_fMIN
        idx_stop=z
        break
      end
    end

    if count_pre==count_f && delay_pre<0 &&count_post>=count_f
      idx_stop=z
      println(count_pre)
      println(count_post)
      break
    end
    if count_post==count_f && delay_pre>=0 &&count_pre>=count_f
      idx_stop=z
      println(count_pre)
      println(count_post)
      break
    end
  end


  return   ww, count_pre, count_post, idx_stop
end
