close all
clear all

% code for simulating the HH equations am using S. Elia paper for reference
% values
tic
%% Some parameters

%time step in s
Nt   =  20000;
dt   =  1e-6;
T    =  0:dt:Nt*dt;
Tinj =  0.002;% time for current injection

%starting voltage and other voltage parameters in V
V0   =  0;
V    =  zeros(1,Nt);
V(1) =  V0;
Vrest= -0.050; % Resting Voltage 
VK   = -0.072;
VNa  =  0.055;
Vl   = -0.050;
%capacitance/area in F/cm^2
Cm   =  1*1e-6;   
%conductivity/area in S/cm^2
G_Na =  0.120; 
G_K  =  0.036;
G_l  =  0.0002; %leak channels
%gating parameters
n0   =  0.32;
h0   =  0.6;
m0   =  0;

n    = zeros(1,Nt);n(1) = n0;
m    = zeros(1,Nt);m(1) = m0;
h    = zeros(1,Nt);h(1) = h0;

% current Injection(A)
I    =  1.2e-4;% from Koch Book

%% An euler implementation
t = 0;
for ii=1:Nt
    t       = t + dt; 
    if(t<Tinj)
          V(ii+1) = V(ii) + dt*1/Cm*(I-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));    
    else
          V(ii+1) = V(ii) + dt*1/Cm*(-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));
    end
    
    %V(ii+1) = V(ii) + dt*1/Cm*(I-(G_K*(n(ii)^4)*(V(ii)-VK) + G_Na*(m(ii)^3)*h(ii)*(V(ii)-VNa) + G_l*(V(ii)-Vl)));     
    V1      = V(ii) - Vrest;% to be used in evaluating gating coefficients
    alpha_n = (0.1-0.01*V1)/(exp(1-0.1*V1)-1);
    beta_n  = 0.125/(exp(0.0125*V1));
    
    alpha_m = (2.5-0.1*V1)/(exp(2.5-0.1*V1)-1);
    beta_m  = 4/(exp(V1/18));
    
    alpha_h = 0.07/(exp(0.05*V1));
    beta_h  = 1/(exp(3-0.1*V1)+1);
    
    n(ii+1) = n(ii) + dt*(alpha_n*(1-n(ii)) - beta_n*n(ii));
    m(ii+1) = m(ii) + dt*(alpha_m*(1-m(ii)) - beta_m*m(ii));
    h(ii+1) = h(ii) + dt*(alpha_h*(1-h(ii)) - beta_h*h(ii));
    
end
    
plot(T,V);

figure
plot(T,m)
hold
plot(T,h)
plot(T,n)
hold

V(1)
V(end)
o =1;
 
toc





     
