function model1rbetak(data,i0,e0)


disp('The decision of r,beta,k using simulated annealing algorithm');

% i0 - the alive infected cases
% e0 - the Potential infected cases
% data - the accmulative cases before using medicine
% r,beta,k - the paras to be decided
% r  0.1 , 0.28  para(1)
% beta 0 ,1  para(2)
% k  1/21 , 1  para(3)


len = length(data)-1;
paras0=[];      % store the best paras
Opt=inf;    % initial the best value 
rand('state',sum(clock));   % set the random seed
for j=1:1000
    para=[(10+18*rand)/100 rand (1+20*rand)/21]; 
    % randperm(n) generates random paras r ,beta , k
    beta2 = 0.75*para(2);
    I = zeros(len,1);
    E = zeros(len,1);
    C = zeros(len,1);
    I(1) = i0;
    E(1) = e0;
    C(1) = i0+e0;
    for i = 1:len,
        beta1 = beta2 + (para(2)-beta2)*exp(-0.015*i);
        E(i+1) = E(i) + beta1*I(i) - para(3)*E(i);
        I(i+1) = para(3)*E(i) + (1-para(1))*I(i);
        C(i+1) = C(i) + para(3)*E(i);
    end
    zz = C(length(C));
    if abs(zz - data(length(C)))>300, % exclude the matlab's error of inf
        para = paras0;
        continue
    end
    cc = data - C;
    ccc = sum(cc.^2);
    if ccc<Opt  % updata the minimal distance
        paras0 = para; % store the paras of the corresponding
        Opt=ccc; 
        D = C;
    end
end
% the above is to generate a better initial paras
e=0.1^35;L=1000000;at=0.999;T=1; 
% set the rse, the maximum steps, Temperature attenuation coefficient, 
% the initial temperature

% the process of simulated annealing
for k=1:L
    % generate new paras by changing a para
    c=ceil(3*rand);
    para = paras0;
    if c==1,
        para(1) = (10+18*rand)/100;
    elseif c==2,
        para(2) = rand;
    else
        para(3) = (1+20*rand)/21;
    end
    % Computational cost function value, thus the sum of square
    beta2 = 0.75*para(2);
    I = zeros(len,1);
    E = zeros(len,1);
    C = zeros(len,1);
    I(1) = i0;
    E(1) = e0;
    C(1) = i0+e0;
    for i = 1:len,
        beta1 = beta2 + (para(2)-beta2)*exp(-0.015*i);
        E(i+1) = E(i) + beta1*I(i) - para(3)*E(i);
        I(i+1) = para(3)*E(i) + (1-para(1))*I(i);
        C(i+1) = C(i) + para(3)*E(i);
    end
    zz = C(length(C));
    if abs(zz - data(length(C)))>300, % exclude the matlab's error of inf
        para = paras0;
        continue
    end
    cc = data - C;
    ccc = sum(cc.^2);
    df = ccc - Opt;
    % rules of acceptance
    if df<0  % it is better than the value we has got 
        paras0 = para;
        Opt=Opt+df;
    elseif exp(-df/T)>rand(1) 
        % accept the current solution by the annealing probability
        paras0 = para;
        Opt=Opt+df;
        D = C;
    end
    T=T*at; % update the temprature
    if T<e  % If it has been almost cooling stop (e ¡æ)
        break;
    end
end

% print the paras and the Opt sum of Square
paras0,Opt

%plot the curve original and in our model

plot(1:length(data),data,'k',1:length(data),D,'r')

    
    

