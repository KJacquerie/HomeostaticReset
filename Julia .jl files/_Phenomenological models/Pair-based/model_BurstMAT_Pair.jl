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


function dx(x::Float64, V::Float64, Vprev::Float64)
  (dt)*-x/tau_p + 1*(Vprev<= 0 && V >0)
end

function dy(y::Float64, V::Float64, Vprev::Float64)
  (dt)*-y/tau_m + 1*(Vprev<= 0 && V >0)
end

function dw(new::Float64, w::Float64, tau::Float64)
  (dt)*(new - w)*1/tau
end



#Sigma_Shouval(Ca::Float64) = a0 - a0.* exp(b1.*(Ca-a1))./(1+exp(b1*(Ca-a1))) + exp(b2.*(Ca-a2))./(1+exp(b2.*(Ca-a2)))



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



#dsSYN(V::Float64, sSYN::Float64) = (dt)*()

function simulateTOY_ncells(ncells::Int64,nEcells::Int64,nIcells::Int64,IappE::Float64,IappI::Float64,TstepEinit::Int64,TstepEfinal::Int64,IstepE::Float64,TstepIinit1::Int64,TstepIinit2::Int64,TstepIinit3::Int64,TstepIinit4::Int64,TstepIfinal::Int64,IstepI1::Float64,IstepI2::Float64,IstepI3::Float64,IstepI4::Float64, gEC::Float64, freq::Float64, delay_pre::Float64, w_init::Float64, SB::Int64)

  # Initial conditions
  V=-60*ones(ncells)
  Vprev=-60*ones(ncells)
  w= w_init
  x::Float64 = 0.
  y::Float64 = 0.
  x_prev::Float64 = 0.
  y_prev::Float64 = 0.
  stop = 0
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
  last = 2
  new = w
  tau = 0.1





  Spkt = zeros(ncells,T)
  lncells::Array{Int64} = ones(ncells,1)

  ww=zeros(convert(Int64,Tdt*dt))
  w = w_init
  idx_w=1

  VV = zeros(Tdt,ncells)

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
      ll = copy(lncells[j])

      # Excitatory cell
      if j==1#(j<=nEcells)
        Iapp = IappE

        z_delayed = z*dt + delay_pre
        if (mod(z_delayed, freq) > freq - duration && z_delayed >=0  && z<TstartE )
          Iappstep = IstepE
        else
          Iappstep = 0.
        end

        if z >= TstartE && z< TstopE
          Iappstep = IstepE2
        end

        V[j] += dV(gNavec_E[j], gKdvec_E[j], glvec_E[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_E[j], gHvec_E[j], gKCavec_E[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_E[j],k1vec_E[j], k2vec_E[j])

        if Vprev[j] < 0. && V[j] >= 0.
          Spkt[j,ll] = t[z]
          lncells[j] += 1
        end


      end # end j = 1


      # Inhibitory cell
      if j==2 #(j > nEcells && j <= nEcells+nIcells)
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

        if Vprev[j] < 0. && V[j] >= 0.
          Spkt[j,ll] = t[z]
          lncells[j] += 1
        end

      end # end j =2

      #  Cortical cell
      if j==3#( j > nIcells+nEcells && j <=ncells)
        Iapp = IappC
        z_delayed = z*dt

        if (mod(z_delayed, freq) > freq - duration && z_delayed >=0 && z<tStartC/dt )
          Iappstep = IstepC
        else
          Iappstep = 0.
        end

        if (z >= tStartC/dt && z< tStopC/dt)
          Iappstep = IstepC2
        end

        V[j] += dV(gNavec_C[j], gKdvec_C[j], glvec_C[j], V[j], mNa[j], hNa[j], mKd[j], mCaT[j], hCaT[j], mH[j], Ca[j], Iapp, Iappstep, gCaTvec_C[j], gHvec_C[j], gKCavec_C[j])
        Ca[j] += dCa(Vprev[j],mCaT[j],hCaT[j],Ca[j],gCaTvec_C[j],k1vec_C[j], k2vec_C[j])

        if Vprev[j] < 0. && V[j] >= 0.
          Spkt[j,ll] = t[z]
          lncells[j] += 1
        end

      end # end j =3

      if j == 1
          V[j] += (dt)*(1/C)*(-gIEGABAA*GABAA[2]*(Vprev[j]+70))
          V[j] += (dt)*(1/C)*(-gIEGABAB*GABAB[2]*(Vprev[j]+85))

      end

      if j == 3
          V[j] += (dt)*(1/C)*(-gEC*AMPA[1]*(Vprev[j]-0))
          V[j] += (dt)*(1/C)*(-gICGABAA*GABAA[2]*(Vprev[j]+70))
          V[j] += (dt)*(1/C)*(-gICGABAB*GABAB[2]*(Vprev[j]+85))
      end


## Plasticity
        if z*dt > TstepEinit

          if (j <=nEcells)
            x_prev = x
            x += dx(x, V[j], Vprev[j])
            if (V[j] >= 0 && Vprev[j]< 0)
              last = 1  #Pré en dernier
            end
          end

          if (j == 3)
            y_prev = y
            y += dy(y, V[j], Vprev[j])
            if (V[j] >= 0 && Vprev[j]< 0)
              last = 0 #Post en dernier

            end
          end


          if last == 0

            if SB == 1
              new = (1-w)*Ap*x_prev + w
            else
              new = Ap*x_prev + w
            end

            w = new

          elseif last == 1

            if SB == 1
              new = -Am*y_prev*w + w
            else
              new = -Am*y_prev + w
            end

            w = new

          elseif last == 2
            w = w
          end
          last = 2
          if (w > wMax)
            w= wMax
          end
          if (w < wMin)
            w= wMin
          end
        else
          w = w
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


    if(mod(z,1/dt)==0)
      ww[idx_w] = copy(w)
      idx_w+=1
    end
    #VV[z,:] = copy(V')
  end
  wf=w

  return  Spkt, ww,wf
end
