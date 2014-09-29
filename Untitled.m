% Stimulation Model

N=100;
dT=0.01;
t=[0:dT:N];

%Currents Injected
%I(1:numel(t))=0;%no current injection (steady state)
%I(1:500)=50; I(501:numel(t))=0; %  curents for the first .5 seconds5and
%zero for the rest
I(1:numel(t))=50; %5uA for the entire simulation time


%Constants
g_K = 36; %mS/cm^2, maximum conductance of potasium
g_Na= 120; %mS/cm^2, maximum conductance of sodium
g_L= 0.3; %mS/cm^2, maximum leakage conductance
E_K= -12; %mV, voltage potential of potasium
E_Na= 115; %mV, voltage potential of Sodium
E_L= 10.6; %mV, voltage leakage potential
Cm=1; %membrane capacitence

%set initial conditions
Vmt=0;
alphaM= (0.1)*((25-Vmt)/(exp((25-Vmt)/10) - 1));
betaM= (4)*exp(-Vmt/18);

alphaN= (.01) * ((10-Vmt)/(exp((10-Vmt)/10) - 1));
betaN= (.125)*exp(-Vmt/80);

alphaH= (.07)*exp(-Vmt./20);
betaH= 1./(exp((30-Vmt)./10)+ 1);

m(1)= alphaM./(alphaM+betaM); % gives the initial value of m at t=0
n(1)= alphaN./(alphaN+betaN); % gives the initial value of n at t= 0
h(1)= alphaH./(alphaH+betaH); % gives the initial value of h at t=0

%loops through each time step 

for tn=1:numel(t)-1
    %calculate the coefficents at each time step
    alphaM(tn)= (0.1).*((25-Vmt(tn))./(exp((25-Vmt(tn))./10) - 1));%calcualtes alphaM at each time step
    betaM(tn)= (4).*exp(-Vmt(tn)/18);%Calcualtes betaM at each time step

    alphaN(tn)= (.01).* ((10-Vmt(tn))/(exp((10-Vmt(tn))./10) - 1)); %calculate alphaM at each time step
    betaN(tn)= (.125).*exp(-Vmt(tn)./80); %calculates betaN at each time step

    alphaH(tn)= (.07).*exp(-Vmt(tn)./20); %calculates alphaH
    betaH(tn)= 1./(exp((30-Vmt(tn))./10)+ 1);% calculates betaH
    %Calculate currents
 I_Na=( m(tn).^3).*g_Na.*h(tn).*(Vmt(tn)-E_Na); %calculates sodium current
 I_K= (n(tn).^4).*g_K.*(Vmt(tn)-E_K) % calculates potassium current
 I_L= g_L.*(Vmt(tn)-E_L) % calculates leakage current
 I_ion= I(tn)-I_K-I_Na- I_L; %calculates the total ionic current
    
    %calculates the derivitives using Euler formula
 Vmt(tn+1) =Vmt(tn)+ dT.*(I_ion/Cm); %calculates the change in the membrane voltage using euler formula
 m(tn+1)= m(tn)+ dT.* (alphaM(tn).*(1-m(tn))-betaM(tn)* m(tn)); %calculates change in m using euler formula
 n(tn+1) = n(tn)+ dT.*(alphaN(tn) .*(1-n(tn))-betaN(tn)* n(tn)); %calculates change in n using euler formula
 h(tn+1)= h(tn) + dT.*(alphaH(tn).*(1-h(tn))-betaH(tn)*h(tn)); % calculates the change in h using eulers formula
 
end 
%plots the membrane potential in mV vs time
Vmt=Vmt-70;
plot(t,Vmt)
title('Membrane Potential')
xlabel('time[ms]')
ylabel('voltage[mv]')
%plots the conductance of sodium and potassium
figure
plot(t,g_K.*(n.^4),t,g_Na.*(m.^3).*h,'r');
title('conductances')
xlabel('time[ms]')
ylabel('Na & K')