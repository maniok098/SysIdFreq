close all
clc
%% 1.1  1.2

nSamples = 100;

sigma_i = 0.5;
sigma_u = 0.5;

i0 = 1;
u0= 1;

% Uk = u0 + sigma_u*randn(nSamples,1);
% Ik = i0 + sigma_i*randn(nSamples,1);

[Uk,Ik] = simulatorVer1(nSamples,u0,i0,sigma_u,sigma_i);

Rk = Uk./Ik;
figure
plot(Rk)

%% 1.2

Ntotal = 3000;
R_SA = zeros(Ntotal,1);
R_LS = zeros(Ntotal,1);
R_EV = zeros(Ntotal,1);
for nSamples = 1:Ntotal

[Uk,Ik] = simulatorVer1(nSamples,u0,i0,sigma_u,sigma_i);

R_SA(nSamples) = 1/nSamples*sum(Uk./Ik);

R_LS(nSamples) = sum(Uk.*Ik)/sum(Ik.^2);

R_EV(nSamples) = sum(Uk)/sum(Ik);


end

% figure
% plot(1:Ntotal,R_SA)
% hold on
% plot(1:Ntotal,R_LS)
% plot(1:Ntotal,R_EV)

figure
plot(1:Ntotal,R_LS)
hold on
plot(1:Ntotal,R_EV)
legend('R_{LS}','R_{EV}')

%% 1.3


%% useful functions

function [Uk,Ik] = simulatorVer1(nSamples,u0,i0,sigma_u,sigma_i)

Uk = u0 + sigma_u*randn(nSamples,1);
Ik = i0 + sigma_i*randn(nSamples,1);


end
