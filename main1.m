close all
clear all

% code for simulating the HH equations am using S. Elia paper for reference
% values. Values according to HH paper. V in mV alpha beta are all in
% msec-1. Capacitance in uF/cm^2. Channel conductance in mS/cm^2. 
tic
%% Some parameters

%time step in ms
dt   =  1e-6;
Ft   =  10;%final time in ms
T    =  0:dt:Ft;
Tinj =  2;% time for current injection
Nt   = length(T);
%starting voltage and other voltage parameters in V
V0   =  -20;
V    =  zeros(1,Nt);
V(1) =  V0;
Vrest= -50; % Resting Voltage 
VK   = -72;
VNa  =  55;
Vl   = -50;
%capacitance/area in uF/cm^2
Cm   =  1;   
%conductivity/area in uS/cm^2
G_Na =  120; 
G_K  =  36;
G_l  =  0.2; %leak channels
%gating parameters
n0   =  0.32;
h0   =  0.6;
m0   =  0;

n    = zeros(1,Nt);n(1) = n0;
m    = zeros(1,Nt);m(1) = m0;
h    = zeros(1,Nt);h(1) = h0;

% current Injection(uA/cm^2)
I    =  10;% from Koch Book
phi  =  1;
flag =  0;

%% An euler implementation
t = 0;
for ii=1:Nt-1
    t       = t + dt;
    %V1      = V(ii) - Vrest;% to be used in evaluating gating coefficients
    if(t<Tinj)
          V(ii+1) = V(ii) + dt*1/Cm*(I-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));
          flag = flag + 1;
    else
          V(ii+1) = V(ii) + dt*1/Cm*(-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));
    end
    
    %V(ii+1) = V(ii) + dt*1/Cm*(I-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));     
    alpha_n = phi*(-0.01*(V(ii)+50))/(exp(-0.1*(V(ii)+50))-1);% as per assignment 2 on IPL
    beta_n  = phi*0.125*(exp(-0.0125*(V(ii)+60)));
    
    alpha_m = phi*(-0.1*(V(ii)+35))/(exp(-0.1*(V(ii)+35))-1);
    beta_m  = 4*phi*(exp(-(V(ii)+60)/18));
    
    alpha_h = phi*0.07*(exp(-0.05*(V(ii)+60)));
    beta_h  = phi/(exp(-0.1*(V(ii)+30))+1);
    
    n(ii+1) = n(ii) + dt*(alpha_n*(1-n(ii)) - beta_n*n(ii));
    m(ii+1) = m(ii) + dt*(alpha_m*(1-m(ii)) - beta_m*m(ii));
    h(ii+1) = h(ii) + dt*(alpha_h*(1-h(ii)) - beta_h*h(ii));
    
end
    
plot(T,V);

figure
plot(T,m,'r')
hold
plot(T,h,'b')
plot(T,n,'y')
hold

V(1)
V(end)
o =1;
 
toc





     
