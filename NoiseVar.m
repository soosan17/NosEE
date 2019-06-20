function sigma_w=NoiseVar(N,P,lambda)
%NoiseVar provides noise variance estimate given size of the noisy
%signal y (N X P) and its corresponding Eigenvalues

%This function is complimentary function for NoSS.m in https://ieeexplore.ieee.org/document/8466044
% N is the length of noisy data
% P is number of measurments / Number of sensors
% lambda are Eigenvalues of data

% Reference Authors:
% Soosan Beheshti (soosan@ee.ryerson.ca) and Saba Sedghizadeh
% CITE:
% https://ieeexplore.ieee.org/document/8466044
% Website: https://www.ee.ryerson.ca/~soosan/
% Code Developement: Saba Sedghizadeh & Younes Sadat-Nejad 
% Contact Info: soosan@ee.ryerson.ca , seyedyouns.sadatneja@ryerson.ca,
% Copy right April 2019
%% ------------------------ obtaining estimated noise variance--------------            
 % Section C of paper the paper "Estimating Uknown Variance"
        if N>P % By maximum gap among the sample variances:
           v_bar_sigma = zeros(1,P);
           v_sigma1 = zeros(1,P);
           v_sigma = zeros(1,P);
           for mi = 1:P
               v_bar_sigma(mi) = sum(lambda(mi:end))/(P-mi+1); %Equation (46)
           end
           for mi = 1:P
                v_sigma1(mi) = (lambda(mi)-v_bar_sigma(mi))^2; 
           end
           for mi = 1:P
               v_sigma(mi) = sum(v_sigma1(mi:end))/(P-mi+1); %Equation (45)
           end
           try
           Vs = v_sigma(1:P-1)./v_sigma(2:P); Vs';
           gVs = find(Vs==max(Vs(1:end))); %Equation (47)
           catch
           Vs = v_sigma(1:P-2)./v_sigma(2:P-1); Vs';
           gVs = find(Vs==max(Vs(1:end-3)));
           end
           sig_w_est = sqrt(mean(lambda(gVs+1:P)));  % Equation (48)
           sigma_w = sig_w_est;
        elseif N <= P
           v_bar_sigma = zeros(1,N);
           v_sigma1 = zeros(1,N);
           v_sigma = zeros(1,N);
           for mi = 1:N
               v_bar_sigma(mi) = sum(lambda(mi:end))/(N-mi+1); %Equation (46) for N data
           end
           for mi = 1:N
                v_sigma1(mi) = (lambda(mi)-v_bar_sigma(mi))^2;
           end
           for mi = 1:N
                v_sigma(mi) = sum(v_sigma1(mi:end))/(N-mi+1);  %Equation (45) for N data
           end
           try
           Vs = v_sigma(1:N-1)./v_sigma(2:N);
           gVs = find(Vs==max(Vs(1:end-1))); %Equation (49)
           catch
           Vs = v_sigma(1:N-2)./v_sigma(2:N-1); Vs';
           gVs = find(Vs==max(Vs(1:end-3)));
           end
           sig_w_est = sqrt(mean(lambda(gVs+1:N)));  % estimated sigma_w
           sigma_w = sig_w_est;
        end