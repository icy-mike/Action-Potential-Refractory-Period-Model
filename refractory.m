clear all
close all
V(1)=0; %initial voltage

time = 25; %testing time in ms
dt = 0.005; %test interval
sizeTime = time/dt; %array size

%constants
C_m = 1;
g_K = 36; g_Na = 120; g_L = 0.3;
V_K = -12; V_Na = 115; V_L = 10.5989;

alpha_n = @(V) [-0.01*V+0.1]/[exp(-0.1*V+1)-1];
beta_n = @(V) 0.125*exp(-0.0125*V);
alpha_m = @(V) [-0.1*V+2.5]/[exp(-0.1*V+2.5)-1];
beta_m = @(V) 4*exp(-V/18);
alpha_h = @(V) 0.07*exp(-0.05*V);
beta_h = @(V) 1/[exp(-0.1*V+3)+1];

n(1) = alpha_n(V)/(alpha_n(V)+beta_n(V)); %initial n m h
m(1) = alpha_m(V)/(alpha_m(V)+beta_m(V));
h(1) = alpha_h(V)/(alpha_h(V)+beta_h(V));

I_m(1:dt:time) = 0;
maxAmp(1:time) = 0;
t = 0:dt:time;

for a = 1:1:time
    I_m(1:40) = 40; %initial pulse
    I_m(40:sizeTime) = 0;
    pulseStart = (a-1)*200+1; %dynamic
    pulseEnd = pulseStart + 39;
    I_m(pulseStart:pulseEnd) = 40;

    for i = 1:sizeTime 
        dndt = @(n) alpha_n(V(i))*(1-n)-beta_n(V(i))*n;
        dmdt = @(m) alpha_m(V(i))*(1-m)-beta_m(V(i))*m;
        dhdt = @(h) alpha_h(V(i))*(1-h)-beta_h(V(i))*h;

        n(i+1) = RK(dndt,n(i),dt);
        m(i+1) = RK(dmdt,m(i),dt);
        h(i+1) = RK(dhdt,h(i),dt);

        INa = g_Na*m(i)^3*h(i)*(V-V_Na);

        dVdt = @(V) [I_m(i)-g_K*n(i)^4*(V-V_K)-g_Na*m(i)^3*h(i)* ...
            (V-V_Na)-g_L*(V-V_L)]/C_m;

        V(i+1) = RK(dVdt,V(i),dt);     
    end
    
    if a == 1
        V_0 = V; %original
    end
    if a > 5
        %difference = V-V_0;
        maxAmp(a) = max(V(1000:end));
    end
    
    subplot(5,5,a)
    plot(t,V)
    xlim([0 25])
    ylim([-25 105])
    
    a
    
end    

figure
subplot(2,1,1)
plot(t,V);
ylabel('Membrane voltage (mV)')
title('Action potential with twin pulse stimulus')
ylim([-25 105])

subplot(2,1,2)
plot(t(1:5000),I_m)
xlabel('Time (ms)')
ylabel('Current stimulus (mA)')
ylim([0 50])

figure
plot(maxAmp)
xlabel('Time (ms)')
ylabel('Membrane voltage (mV)')
title('Second AP amplitude')
