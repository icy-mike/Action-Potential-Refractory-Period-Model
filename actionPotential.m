clear all
close all
V(1)=0; %initial voltage

time = 50; %testing time in ms
dt = 0.05; %test interval
sizeTime = time/dt; %array size

%constants
I_m = 14; %current stimulation, uA/cm^2
C_m = 1; %membrane capacitance, uF/cm^2
g_K = 36; g_Na = 120; g_L = 0.3; %max ionic leakage conductance, mmho/cm^2
V_K = -12; V_Na = 115; V_L = 10.5989; %Nernst potentials (offset from resting voltage), mV

alpha_n = @(V) [-0.01*V+0.1]/[exp(-0.1*V+1)-1];
beta_n = @(V) 0.125*exp(-0.0125*V);
alpha_m = @(V) [-0.1*V+2.5]/[exp(-0.1*V+2.5)-1];
beta_m = @(V) 4*exp(-V/18);
alpha_h = @(V) 0.07*exp(-0.05*V);
beta_h = @(V) 1/[exp(-0.1*V+3)+1];

n(1) = alpha_n(V)/(alpha_n(V)+beta_n(V)); %initial n m h
m(1) = alpha_m(V)/(alpha_m(V)+beta_m(V));
h(1) = alpha_h(V)/(alpha_h(V)+beta_h(V));

for i = 1:sizeTime 
    dndt = @(n) alpha_n(V(i))*(1-n)-beta_n(V(i))*n;
    dmdt = @(m) alpha_m(V(i))*(1-m)-beta_m(V(i))*m;
    dhdt = @(h) alpha_h(V(i))*(1-h)-beta_h(V(i))*h;
 
    n(i+1) = RK(dndt,n(i),dt);
    m(i+1) = RK(dmdt,m(i),dt);
    h(i+1) = RK(dhdt,h(i),dt);
    
    GNa(i) = g_Na*m(i)^3*h(i);
    GK(i) = g_K*n(i)^4;

    dVdt = @(V) [I_m-g_K*n(i)^4*(V-V_K)-g_Na*m(i)^3*h(i)* ...
        (V-V_Na)-g_L*(V-V_L)]/C_m;
    
    V(i+1) = RK(dVdt,V(i),dt);     
end
t = 0:dt:time;
t1 = t(100:500);
figure(1)
plot(t,-70+V);
% hold on
% plot(t(1:5000),INa*50);
xlabel('Time (ms)')
ylabel('Membrane voltage (mV)')
title('Action Potential')

figure(2)
plot(t,n)
hold on
plot(t,m,'g')
plot(t,h,'r')
xlabel('Time (ms)')
ylabel('Constant values')
title('Conductance constants')

figure(3)
plot(t(1:1000),GK)
hold on 
plot(t(1:1000),GNa,'g')
xlabel('Time (ms)')
ylabel('Conductance (mmho)')
title('Na^+ and K^+ conductance')

figure(4)
subplot(3,1,1)
plot(t1,V(100:500));
ylabel('Membrane voltage (mV)')
title('Action Potential')
xlim([5 25])

subplot(3,1,2)
plot(t1,n(100:500))
hold on
plot(t1,m(100:500),'g')
plot(t1,h(100:500),'r')
ylabel('Constant values')
title('Conductance constants')
xlim([5 25])

subplot(3,1,3)
plot(t(100:500),GK(100:500))
hold on 
plot(t(100:500),GNa(100:500),'g')
xlabel('Time (ms)')
ylabel('Conductance (mmho)')
title('Na^+ and K^+ conductance')
xlim([5 25])
