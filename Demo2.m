%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Noss is Mean Squared Eigenvalue Error method used for obtaining optimum
% Number of Source Signals (NoSS): https://ieeexplore.ieee.org/document/8466044
% Reference Authors:
% Soosan Beheshti (soosan@ee.ryerson.ca) and Saba Sedghizadeh
% CITE:
% "Number of Source Signal Estimation by the Mean Squared Eigenvalue Error."
% IEEE Transactions on Signal Processing 66.21 (2018): 5694-5704.
% https://ieeexplore.ieee.org/document/8466044
% Website: https://www.ee.ryerson.ca/~soosan/
% Code Developement: Saba Sedghizadeh & Younes Sadat-Nejad 
% Contact Info: soosan@ee.ryerson.ca , seyedyouns.sadatneja@ryerson.ca,
% Copy right April 2019


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Demo of detecting Number of Source Signals (NoSS) from their linear mixture
clear , clc ,close all
P = 25;            % number of observations or sensors
q = 20;            % number of parameters to be detected 
A = zeros(P,q);    % The matrix created as in Dr. Beheshti's 2010 paper
LL = randn(1,q+P);
for i=1:P
    A(i,:) = LL(i:i+q-1);
end
noise_var = 200;
N= 90;             %Data Length
s = zeros(q,N); y = zeros(P,N); y_bar = zeros(P,N);
Mm = 15;
s(1:q-Mm,1:N) = randn(q-Mm,N);  % clean parameters, the remaining 10 
                                % parameters are set to zero 
                                %(the parameters are not constant during the time)
                                % Here the model order is (q-Mm) or num of sources
sigma_w = sqrt(noise_var);  % squared roots of noise variance
    
%------------------ obtaining the noisy observations----------------------
            w = sigma_w*randn(P,N);         % noise matrix   
            W = (1/N)*(w)*(w)';             % noise covariance matrix
            y_bar(:,1:N) = A*s(:,1:N);      % noiseless observations
            y = y_bar + w;                  % noisy observations
%------------------------------ PLOTS ------------------------------------
subplot(2,1,1)
plot(1:N,y_bar);
title('Source Signals (Uknown)')
ylabel('Amplitude')
xlabel('Number of Samples')
subplot(2,1,2)
plot(1:N,y);
title('Noisy Mixutre of Source Signals')
ylabel('Amplitude')
xlabel('Number of Samples')


%% Using NoSEE algorithm to detect number of sources (Expected solution is 3)
m=NoSEE(y); %  SAME AS m=NoSS(y,5,5,1000);
disp('** NoSS Result **'); 
disp(m)
disp('** True Result **'); 
disp(q)

% Get Number Of Source Signals AND SVD in MATLAB 
      %[m,U,S,V]=NoSEE(Y);
     
