function w_predict = predictw(demo,expm)

tau_w    = 520.76129e3;
gamma_p  = 597.08922;
gamma_d  = 137.7586;
theta_p = 2.009289;
theta_d = 1.0;


Omega_p = gamma_p/(gamma_p+gamma_d);
Omega_d = 0;
Omega_0 = 0;

tauw_p = tau_w /(gamma_p + gamma_d);
tauw_d = tau_w / gamma_d;
tauw_0 = 0;

zeta = tauw_d / tauw_p; 

w_predict = zeros(1,size(demo,1));

if expm == "Control"
    for i=1:1:size(demo,1)
        tot = demo(i,4) - demo(i,1); % if we don't consider when < theta_d (total - r) 
        alpha_d = demo(i,2)/tot;
        alpha_p = demo(i,3)/tot;

        w_predict(i) = Omega_p*((alpha_p*zeta)/(alpha_p*zeta+alpha_d));    
    end
end

if expm == "DD_rel_half"
    r=0.5; %1 if rel, 0.5 if rel_half
    for i=1:1:size(demo,1)
        tot = demo(i,4);
        alpha_d = demo(i,2)/tot;
        alpha_p = demo(i,3)/tot;
        w_predict(i) = alpha_p*(r*((gamma_p-gamma_d)/tau_w))+alpha_d*(-r*gamma_d/tau_w);
    end
end
    
if expm == "DD_rel"
    r=1;
    for i=1:1:size(demo,1)
        tot = demo(i,4);
        alpha_d = demo(i,2)/tot;
        alpha_p = demo(i,3)/tot;
        w_predict(i) = alpha_p*(r*((gamma_p)/tau_w))+alpha_d*(-r*gamma_d/tau_w);
    end
end 

if expm == "Control/checkFD"
    for i=1:1:size(demo,1)
        tot_stop = demo(i,4) - demo(i,1); % if we don't consider when < theta_d (total - r) 
        alpha_d = demo(i,2)/tot;
        alpha_p = demo(i,3)/tot;
        w_predict(i) = Omega_p*((alpha_p*zeta)/(alpha_p*zeta+alpha_d));    
    end
end

end