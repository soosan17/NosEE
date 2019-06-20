function [m,U,S,V]=NoSEE(y,alpha,beta,itt)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% NosEE is Number of Source estimate by Eigenvalue Error
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
% ---- INPUT ---
%  y : Noisy Measurments : y = y_bar + noise
%  alpha , beta : Probabilistic values (OPTIONAL : Both equal 4 as default)
%  itt : Number of iterations for Monte Carlo Loop ( OPTIONAL: 1000 as default)
%     Higher value for itt provides better estimation of noise variance
%     Lower  value provides faster result (Lower accuracy)
% ---- OUTPUT ---
% m: Optimum Number of Source Signals (NoSS) for the noisy signal y
% U,S,V]: result of SVD from MATLAB svd Algorithm 
% ---------------
%   Class Support
%   The input signal ( y ) is P x N matrix , consist of real numbers where:
% P is number of measured signals (sensors)
% N is length of signal
%   alpha and beta are float numbers
%   itt is real integer
% ---------------
%   Example of Use:
%   [m,U,S,V]=NoSEE(y) : Provides Number of Sources in y and U,S,V of
%   MATLAB "svd" function
%   m=NoSEE(y)         : Provides Number of Sources in y 
%                        alpha=beta=4 ,itt = 1000 by default
%   m=NoSEE(y,alpha,beta,itt): Privies Number of Sources in y using
%   probalistic values of alpha,beta and numbter of itteration "itt"
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%



 validateattributes(y,{'numeric'},...
    {'real','finite','nonsparse'},2)

if nargin < 2
    alpha=5;
    beta=5;
    itt=1000;
end
if ~isempty(y)

[P,N]=size(y);
%P Number of observations (Number of Sensors)
%N Lenght of observation
                    % run itteration for tt times
%----- obtaining the noisy sample covariance matrix and its eigenvalues---
            C_hat = (1/N)*((y)*(y)');   %Equation (3)   % sample covariance matrix
            lambda = svd(C_hat);      %Eigenvalues of C_hat

if nargout == 4 
            [U,S,V] = svd(y);
end
for tt=1:itt   
%----- obtaining the noisy sample covariance matrix and its eigenvalues---
sigma_w=NoiseVar(N,P,lambda); % Section C of paper the paper "Estimating Uknown Variance"
% ------------------------- Initialization for MSEE ---------------------
            svCov_y = lambda;
            svCov_y_m = zeros(size(svCov_y));
            Xm = zeros(size(svCov_y));
            Delta_m = zeros(size(svCov_y));
            U_Zm = zeros(size(svCov_y));
            L_Zm = zeros(size(svCov_y));
            Sw = zeros(size(svCov_y));
% ------------------------- Main loop of MSEE ---------------------
      for m = 1:P
            Sw = zeros(size(svCov_y));
            VarS(m) = 2*((sigma_w)^4)/N + 4*((sigma_w)^2)*mean(var(y))/N;
            Sw(1+m:P) = (sigma_w)^2;
            svCov_y_m(1:m) = svCov_y(1:m);
            Xm(m) = norm(svCov_y - svCov_y_m - Sw)^2;  % available eigenvector error
            L_alpha(m) = N*((P-m)/N - Xm(m)/(sigma_w)^4)/sqrt(2*(P-m));
 % ************** Estimated Delta and its Upper and Lower bounds **********           
            mw(m) = (P-m)*VarS(m); %Equation (29)
            Delta_mest(m) = Xm(m) - mw(m);     %Equation          % available Delta_m
            Km(m) = 2*alpha*sqrt(mw(m)*((alpha^2/(P-m)-(1/2))*mw(m) + Xm(m))/(P-m)); %Equation (31)
            U_Delta_mest(m) = Xm(m) + (2*alpha^2/(P-m) - 1)*mw(m) + Km(m); %Equation (27)
            L_Delta_mest(m) = Xm(m) + (2*alpha^2/(P-m) - 1)*mw(m) - Km(m); %Equation (28)
                     
% ************* Order estimate from Upper bound of Zm  *******************
            EZm(m)   = m*VarS(m) + m*sigma_w^4 + U_Delta_mest(m); %Equation (22)
            VarZm(m) = 2*m*VarS(m)*(VarS(m) + 2*(sigma_w)^4); %Equation(23) 
% *************** Upper bound and Lower bound for Zm  ********************
            U_Zm(m) = m*VarS(m) + m*sigma_w^4 + U_Delta_mest(m) + beta*sqrt(VarZm(m)); %Equation (34)
      end
    min_MSEE(tt) =find(U_Zm==min(U_Zm));
    
     end
    m = mean(min_MSEE);

end
