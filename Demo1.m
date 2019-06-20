%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%Noss is Mean Squared Eigenvalue Error method used for obtaining optimum
%Number of Source Signals (NoSS): https://ieeexplore.ieee.org/document/8466044
% Reference Authors:
% Soosan Beheshti (soosan@ee.ryerson.ca) and Saba Sedghizadeh
% CITE:
%"Number of Source Signal Estimation by the Mean Squared Eigenvalue Error."
%IEEE Transactions on Signal Processing 66.21 (2018): 5694-5704.
% https://ieeexplore.ieee.org/document/8466044
% Website: https://www.ee.ryerson.ca/~soosan/
% Code Developement: Saba Sedghizadeh & Younes Sadat-Nejad 
% Contact Info: soosan@ee.ryerson.ca , seyedyouns.sadatneja@ryerson.ca,
% Copy right April 2019
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Demo of detecting Number of Source Signals (NoSS) from their linear mixture
clear, clc , close all
N   = 5000;             % Number samples
NoisePower = 0.1;         % Power of additive Noise
d   = 7;                % Number of Observations (After Linear Mixture)
r   = 3;                % Number Source Signals (Before Linear Mixture-Uknown)
% Hence
% Generating y_bar (Unknown Signals)
t = @(nn,P) linspace(0,1,nn) * 2 * pi * P;
y_bar(1,:) = sin(t(N,2));            % Sine Signal
y_bar(2,:) = cos(t(N,8));            % Cosine Signal
y_bar(3,:) = sign(sin(t(N,5)));      % Square Wave Signal
% Adding noise to the signals before mixing 
sigma= NoisePower; %Using SNR value given previously
y_noisy = y_bar + sigma * randn(size(y_bar));

% Generating the Mixed Signals
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
A = normRows(rand(d,3));
Y = A * y_noisy; 

% Plots to illustrate mixing source signals and signals after mixing:
figure;
subplot(3,1,1)
plot(1:N,y_bar)
xlabel('Sample Points')
ylabel('Amplitude')
title('Source Signals (3 Signals)')
legend('Sine','Cosine','Square')
subplot(3,1,2)
plot(1:N,y_noisy)
xlabel('Sample Points')
ylabel('Amplitude')
title('Noisy Source Signals (3 Signals)')
legend('Sine','Cosine','Square')
subplot(3,1,3)
plot(1:N,Y)
xlabel('Sample Points')
ylabel('Amplitude')
title('Noisy Signals After Mixing (7 signals)')

%% Using NoSEE algorithm to detect number of sources (Expected solution is 3)
m=NoSEE(Y); %  SAME AS m=NoSS(y,5,5,1000);
disp('** The Number of Source Signals with 1000 itteration and alpha,beta=4 is: **'); 
disp(' ')
disp(m)

% Get Number Of Source Signals AND SVD in MATLAB 
      [m,U,S,V]=NoSEE(Y);
