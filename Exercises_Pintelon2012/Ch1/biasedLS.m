close all
clc

% noise,  zero mean
sigma_i = 0.5;
sigma_u = 0.5;

% real parameter: R = u0/i0 = 1
i0 = 1;
u0= 1;

% Uk = u0 + sigma_u*randn(nSamples,1);
% Ik = i0 + sigma_i*randn(nSamples,1);

Ntotal = 4000;

R_LS = zeros(Ntotal,1);  % Least Square
R_EV = zeros(Ntotal,1);  % errors-in-variables
for nSamples = 1:Ntotal

[Uk,Ik] = simulatorVer1(nSamples,u0,i0,sigma_u,sigma_i);

R_LS(nSamples) = sum(Uk.*Ik)/sum(Ik.^2);  % derive it :-)

R_EV(nSamples) = sum(Uk)/sum(Ik);


end

figure
plot(1:Ntotal,R_LS)
hold on
plot(1:Ntotal,R_EV)
legend('R_{LS}','R_{EV}')

%% useful functions

function [Uk,Ik] = simulatorVer1(nSamples,u0,i0,sigma_u,sigma_i)

Uk = u0 + sigma_u*randn(nSamples,1);
Ik = i0 + sigma_i*randn(nSamples,1);


end
