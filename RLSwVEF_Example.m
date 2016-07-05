%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
% this programe is used to estimate the parameter of 
%                     0.5
% G(s) = -----------------    
%               s^2 + s + 1 
%                0.009335 z + 0.008732
% G(z) = ---------------------------------------------           ; T = 0.2
%                z^2 - 1.783 z + 0.8187
%% Plot option
% if Plot = 0  --> no plot needed
% if Plot = 1  -->  plot needed
Plot=1;
%% input signal
T=0.2;
t=0:T:50;
u=2*exp(-0.01*t).*sin(0.1*t);
%% output signal
y(1)=0;
for k=1:length(u)-1
    [ y_output ] = OutputEstimation( [1 -1.783  0.8187], [0.009335  0.008732], 1, u(1:k+1), y(1:k), k+1 );
    y(k+1)=y_output;
end
%% 1 st order estimation
[ Gz1 ] = RecursiveLeastSquarewithVaryingExponentialForgetting ( u, y, 1, 1, 0,  T, [Plot 1 2 3] )
%% 2 st order estimation
[ Gz2 ] = RecursiveLeastSquarewithVaryingExponentialForgetting ( u, y, 2, 1, 1, T, [Plot 4 5 6] )
%% 3 st order estimation
[ Gz3 ] = RecursiveLeastSquarewithVaryingExponentialForgetting ( u, y, 3, 1, 1, T, [Plot 7 8 9] )
