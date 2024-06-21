% Gillespie algorithm for the cultural shift model
% simulates the gillespie algorithm n times

function [] = cutural_shift_gillespie



n=1; % number of realisations

for j=1:n % n simulations at a parameter value


 clearvars -except j n k time_of_emerg number_of_mutations p_thresh prob_emerg prob_emerg_low timestep time_of_first_mut M_new M_dead


% N1 and N2 parameters

  rho = 0.005;          % fitness advantage from culture shift [-0.6 0 1]
  alpha1 = 3*10^(-5);   % Soc-learning N1 -> N2 [5*10^(-6) to 5*10^(-5)]
  alpha2 = 0.01;        % Abandonment of practice N2,N2->N1 [0.05 0.03 0.1]
  alpha3 = 0*alpha1;  % Soc-learning N2 -> N3
  alpha4 = 0*alpha2;    % Abandonment N3 -> N1 (general model in SI)
  alpha5 = 0*alpha1;    %Soc-leaning N1 -> N3 (general model in SI)
  alpha6 = 0*alpha2;    %abandonment N3 -> N2 (general model in SI)
  gamma = 5*10^(-5);    % mutation from N2 to N3 5*10^(-5)
  gamma_N2 = 0*10^(-4); % mutation from N1 to N2 (we do not use this in the model)
  lambda = 1;           % relative fitness benefit for N3
  sig = 0.1;            % relative risk of practice from N3
  a = 0;                % hill function coefficient (nonlinear microbe interaction in SI)
  K=10^5;               % half saturation growth rate (nonlinear microbe interaction in SI)
  D = 0;                % conformity bias coefficient
  death = 0;            % death from infection (used in SI version)
 
  rm = 10^(-5);         % M growth 
  dm = 0.1;             % M death rate 
  rp = 1.5*10^(-5);     % P growth
  sp = 10^(-5);         % removal of P
  u = 5*10^(-3);        % mutation M to P 
 
  

% Setup dynamical array and initial condition
  t(1) = 0;
  N1(1) = 9999;
  N2(1) = 1;
  N3(1) = 0;
  M(1) = 0;
  P(1) = 0;
  mut(1) = 0;
  M_birth(1) = 0;
  M_death(1) = 0;

 

tend = 500; % End time of simulation
i = 1;

while t(end)<tend


% Rates

    n2 = N2(i)/(N1(i)+N2(i));
    if isnan(n2)
        n2=1;
    end


    rates = zeros(1,19);
    
    rates(1) = sp*P(i)*(N1(i)+N2(i)+N3(i));                 % P death
    rates(2) = alpha1*N1(i)*N2(i)*(1+D*(1-n2)*(2*n2-1));    % Social learning N1 -> N2
    rates(3) = rho*N2(i);                                   % Reproductive advantage for N2
    rates(4) = alpha4*N3(i);                                % Abandonment N3 -> N1
    rates(5) = rm*N2(i)*(N2(i)^a/(N2(i)^a+K^a-1))...
        + rm*N3(i)*sig*(N3(i)^a/(N3(i)^a+K^a-1));           % M growth
    rates(6) = dm*M(i);                                     % M death
    rates(7) = u*M(i);                                      % Mutation M -> P
    rates(8) =  rp*(N1(i)+N2(i)+N3(i))*P(i);                % P growth
    rates(9) = alpha2*N2(i);                                % Abandonment of practice N2 -> N1
    rates(10) = rho*lambda*N3(i);                           % Reproductive advantage for N3
    rates(11) = alpha3*N2(i)*N3(i)*(1+D*(1-n2)*(2*n2-1));   % Social learning N2 -> N3
    rates(12) = gamma*N2(i);                                % discovery of N3
    rates(13) = alpha5*N1(i)*N3(i);                         % Social learning N1 -> N3
    rates(14) = death*N2(i)*P(i);                           % Death due to infection N2
    rates(15) = gamma_N2*N1(i);                             % mutation from N1 to N2 (we do not use this in the model)
    rates(16) = alpha6*N3(i);                               % Abandonment N3 -> N2
    rates(17) = death*N1(i)*P(i);                           % mutation from N1 to N2 (we do not use this in the model)
    rates(18) = death*N3(i)*P(i);                           % mutation from N1 to N2 (we do not use this in the model)
    rates(19) = 0*dm*P(i);                                  % death rate of P (used in SI version)
    rate_sum = sum(rates);

  % conditions to end simulations and end muation

    if rate_sum==0
        break;
    end

    if P(i)==100
            break
    end

%     if N3(i)>=1 % we can use this if we want a only single individual to invent the alt practice
%         gamma = 0;
%     end

%     if N1(i)<0 || N2(i)<0 || N3(i)<0
%         break
%     end


    r1 = rand(1);
    tau = (1/rate_sum)*log(1/r1);   %increment of time 

    t(i+1) = t(i) + tau;            % next time point when an event occurs

    r2 = rand(1);                   % next event 
    

%rates 1
    if r2*rate_sum <= rates(1) %&& r2*rate_sum <= rates(1) + rates(2)
         N1(i+1) = N1(i);
         N2(i+1) = N2(i)-1;
         N3(i+1) = N3(i);
         M(i+1) = M(i);
         P(i+1) = P(i) - 1;
         mut(i+1) = mut(i);
         M_birth(i+1) = M_birth(i);
         M_death(i+1) = M_death(i);
%rates 2
    elseif r2*rate_sum > rates(1)  && r2*rate_sum <= rates(1) + rates(2) 
        N1(i+1) = N1(i) - 1;
        N2(i+1) = N2(i) + 1;
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 3
    elseif r2*rate_sum > rates(1) + rates(2) && r2*rate_sum <= rates(1) + rates(2) + rates(3) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i) + 1;
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 4
    elseif r2*rate_sum > sum(rates(1:3))  && r2*rate_sum <= sum(rates(1:4)) 
        N1(i+1) = N1(i) + 1;
        N2(i+1) = N2(i);
        N3(i+1) = N3(i) - 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 5
    elseif r2*rate_sum > sum(rates(1:4))  && r2*rate_sum <= sum(rates(1:5)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i); 
        M(i+1) = M(i) + 1;
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i)+1;
        M_death(i+1) = M_death(i);

%rates 6
    elseif r2*rate_sum > sum(rates(1:5)) && r2*rate_sum <= sum(rates(1:6)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i);
        M(i+1) = M(i) - 1;
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i)-1;

%rates 7
    elseif r2*rate_sum > sum(rates(1:6)) && r2*rate_sum <= sum(rates(1:7))
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i); 
        M(i+1) = M(i) - 1;
        P(i+1) = P(i) + 1;
        mut(i+1) = mut(i) + 1;
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 8
    elseif r2*rate_sum > sum(rates(1:7))  && r2*rate_sum <= sum(rates(1:8)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i); 
        M(i+1) = M(i);
        P(i+1) = P(i) + 1;
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 9
    elseif r2*rate_sum > sum(rates(1:8)) && r2*rate_sum <= sum(rates(1:9)) 
        N1(i+1) = N1(i) + 1;
        N2(i+1) = N2(i) - 1;
        N3(i+1) = N3(i); 
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);
    
%rates 10
   elseif r2*rate_sum > sum(rates(1:9))  && r2*rate_sum <= sum(rates(1:10))
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i) + 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);
    
    
%rates 11
    elseif r2*rate_sum > sum(rates(1:10)) && r2*rate_sum <= sum(rates(1:11))
        N1(i+1) = N1(i);
        N2(i+1) = N2(i) - 1;
        N3(i+1) = N3(i) + 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 12
    elseif r2*rate_sum > sum(rates(1:11)) && r2*rate_sum <= sum(rates(1:12)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i) - 1;
        N3(i+1) = N3(i) + 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 13
    elseif r2*rate_sum > sum(rates(1:12)) && r2*rate_sum <= sum(rates(1:13)) 
        N1(i+1) = N1(i) - 1;
        N2(i+1) = N2(i);
        N3(i+1) = N3(i) + 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 14
    elseif r2*rate_sum > sum(rates(1:13)) && r2*rate_sum <= sum(rates(1:14)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i) - 1;
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

%rates 15
    elseif r2*rate_sum > sum(rates(1:14)) && r2*rate_sum <= sum(rates(1:15)) 
        N1(i+1) = N1(i) - 1;
        N2(i+1) = N2(i) + 1;
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);
 
  %rates 16
    elseif r2*rate_sum > sum(rates(1:15)) && r2*rate_sum <= sum(rates(1:16)) 
        N1(i+1) = N1(i);
        N2(i+1) = N2(i) + 1;
        N3(i+1) = N3(i) - 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);  

   %rates 17
    elseif r2*rate_sum > sum(rates(1:16)) && r2*rate_sum <= sum(rates(1:17))
        N1(i+1) = N1(i) - 1;
        N2(i+1) = N2(i);
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

    %rates 18
    elseif r2*rate_sum > sum(rates(1:17)) && r2*rate_sum <= sum(rates(1:18))
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i) - 1;
        M(i+1) = M(i);
        P(i+1) = P(i);
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);

     %rates 19
     elseif r2*rate_sum > sum(rates(1:18)) && r2*rate_sum <= sum(rates(1:19))
        N1(i+1) = N1(i);
        N2(i+1) = N2(i);
        N3(i+1) = N3(i);
        M(i+1) = M(i);
        P(i+1) = P(i) - 1;
        mut(i+1) = mut(i);
        M_birth(i+1) = M_birth(i);
        M_death(i+1) = M_death(i);   
    end   
    
i = i+1;

end


ind = find(P>=100,1,'first'); % Threshold for emergence

if isempty(ind)

    time_of_emerg(1,j) = NaN;
    number_of_mutations(1,j) = NaN;
else 
   
    time_of_emerg(1,j) = t(ind);

    % number of muataions at emergence
    number_of_mutations(1,j) = mut(ind);
end

M_new(j) = M_birth(end);
M_dead(j) = M_death(end);

ind_mut = find(mut==1,1,'first');

if isempty(ind_mut)
   time_of_first_mut(1,j)=NaN;
else
   time_of_first_mut(1,j) = t(ind_mut);
end
time_of_first_mut;

mutations = [0 diff(mut)]; % create and array with the mutations at each timepoint


%figures where all populations are represented in subplots
figure

% subplot(4,1,1)
% plot(t,mutations,'LineWidth',3)
% set(gca,'FontSize',15,'fontweight','bold')
% title('Mutations')
% size(N1)
% size(N2)
% size(N3)
% size(t)
subplot(2,1,1)
plot(t,N1,t,N2,t,N3,'LineWidth',2)
%plot(t,N2,'LineWidth',3)
%set(gca,'FontSize',20,'fontweight','bold')
ylabel('Human populations')
legend('N_1','N_2','N3','Location','eastoutside')
%set(legend, 'Box', 'off');
xlim([0 t(end)])
box off

subplot(2,1,2)
semilogy(t,M,'k',t,P,'m','LineWidth',2)
%set(gca,'FontSize',20,'fontweight','bold')
ylabel('Microbe/Pathogen population')
xlabel('Time (years)')
legend('M','P','Location','eastoutside')
%set(legend, 'Box', 'off');
xlim([0 t(end)])
ylim([0.95 100])
box off


end
end
