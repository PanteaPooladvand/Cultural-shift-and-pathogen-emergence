% Runs the cultural_shift_ODE model. This solver can be used to run the
% model at the set paramter values.
% The solver can also be used to do a sensitivity analysis using the arrays
% labelled as pp.


function [time,pop] = solver_cultural_shift_ODE

time = {};
pop = {};

%pp = 1; % use if running the model at the set paramter values. Remove pp(i) from rho or alpha_3 and replace with the given parameter values
pp = [0.005 0.01 0.02 0.04]; % use for rho sensitivity analysis
%pp = [0 0.5 1 2]; % use alpha3/alpha1 sensitivity analysis
%pp = [0 1 10 100];

tspan = [0 1000];

param_names = {'\alpha_3/\alpha_1=' '\rho =' }; % parameters in legend
name = 2;

for i = 1:length(pp)
i


% N1 and N2 parameters

  r1 = 0;       % N1 growth
  rho = pp(i);%0.005;  % fitness advantage of practice
  alpha1 = 3*10^(-5); % Soc-learning N1->N2
  
  % M and P parameters

  rm = 10^(-5);        % M growth
  dm = 0;            % M death rate 
  rp = 1.5*10^(-5);    % P growth 
  sp = 10^(-5);        % fitteness due to P
  u = 5*10^(-3);       % mutation from M to P 
  alpha2 = 0.1;       % abandonment N2 -> N1
  alpha4 = 0*alpha2;   % abandonment N3 -> N1
  alpha3 = 1/10*alpha1; % social learning N2 -> N3 set this at 0.1*alpha1 when rho is preturbed
  alpha5 = 0*alpha1;   % socal learnin N1 -> N3
  alpha6 = 0*alpha2;   % abandonment N3 -> N2
  sigma = 0.1;           % relative risk of alternative practice N3 
  lambda = 1;          % relative reprod advantage
  D = 0;
  gamma = 5*10^(-5);   % discovery of N3
  a = 0;
  K = 10^5;
  delta = 0;
  delta1 = delta;
  delta2 = delta;
  delta3 = delta;

param1 = [r1 rho alpha1];
param2 = [rm dm rp sp u alpha2 alpha4 alpha3 lambda alpha5 alpha6 sigma a K D gamma delta1 delta2 delta3];

y0 = [9999 1 0.9 0.9 0];

options = odeset('RelTol',1.e-6);

sol = ode45(@(t,y)cultural_shift_ODE(t,y,param1,param2),tspan, y0,options);


t = sol.x;
y = sol.y;

time{i} = t;
pop{i} = y;
max(y(4,:));
end

figure
hold on

% subplot(1,2,2)
for i = 1:length(pp)
    
    N1 = pop{i}(1,:); N2 = pop{i}(2,:); N3 = pop{i}(5,:); % extract the N populations
    sty = {'-',':','-.','--'};
    txt = [sprintf('%s',param_names{name},num2str(pp(i)))]; % set legend
    plot(time{i},N2 ,string(sty(i)),'DisplayName',txt,'color',[0.8500 0.3250 0.0980],'Linewidth',2) % plot populations
    set(gca,'FontSize',20,'fontweight','bold')
    box off
    axis square
    ylim([0 6.7*10^4])
    xlim([0 300])
    xlabel('Time (years)')
    ylabel('Frequency in original practice')
end

legend show
hold off

%% Use this figure for plotting populations with only a single loop
N1 = pop{1}(1,:); N2 = pop{1}(2,:); N3 = pop{1}(5,:); % extract the N populations
M = pop{1}(3,:); P = pop{1}(4,:);

figure

subplot(2,1,1)
plot(t,N1,t,N2,t,N3,'LineWidth',2)

ylabel('Human populations')
legend('N_1','N_2','N3','Location','eastoutside')
set(legend, 'Box', 'off');
xlim([0 t(end)])
box off

subplot(2,1,2)
semilogy(t,M,'k',t,P,'m','LineWidth',2)

ylabel('Microbe/Pathogen population')
xlabel('Time (years)')
legend('M','P','Location','eastoutside')
set(legend, 'Box', 'off');
xlim([0 t(end)])
ylim([0.95 100])
box off




