% This is the general ODE cultural shift model.
% It includes all populations N1 to N3 and interactions between them.
% Use the solver_cultural_shift_ODE to change paramter values and to run the model.

function dydt = cultural_shift_ODE(t,y,param1,param2)

  N1 = y(1);
  N2 = y(2);
  M = y(3);
  P = y(4);
  N3 = y(5);

  % Parameters

  r1 = param1(1);            % N1 growth
  rho = param1(2);           % N2 growth
  alpha1 = param1(3);        % Soc-learning N1 -> N2

  rm = param2(1);            % M growth 
  dm = param2(2);            % M death rate 
  rp = param2(3);           % P growth from N1
  sp = param2(4);         % fitteness due to P
  u = param2(5);            % mutation N to P

  alpha2 = param2(6);        % Abandonment N2 -> N1
  alpha4 = param2(7);        % Abandonment N3 -> N1
  alpha3 = param2(8);        % Soc-learning N1 -> N2
  lambda = param2(9);        % relative fitness advantage N3
  alpha5 = param2(10);       % Soc-learning N1 -> N3
  alpha6 = param2(11);       % Abandonment N3 -> N2
  sigma = param2(12);        % relative risk on practice from N3 (microbe M growth)
  a = param2(13);            % hill coefficient
  K = param2(14);            % half-saturation constant
  D = param2(15);            % frequency dependance coefficient for conformity
  gamma = param2(16);        % discovery of N3 practice
  delta1 = param2(17);       % death of N1 due to pathogen P
  delta2 = param2(18);       % death of N2 due to pathogen P
  delta3 = param2(19);       % death of N2 due to pathogen P
  n2 = N2./(N1+N2);          % if using conformity bias 


  % ODEs
  
  dydt(1) = r1*N1  - alpha1*N1*N2 - alpha5*N1*N3 + alpha2*N2 + alpha4*N3 - delta1*sp*N1*P;                      % N1 population
  dydt(2) = N2*rho + alpha1*N1*N2 - alpha2*N2 - alpha3*N2*N3 + alpha6*N3 - delta2*sp*N2*P;                      % N2 population
  dydt(3) = rm*N2*(N2^a/(N2^a+K^a-1)) + sigma*rm*N3*(N3^a/(N3^a+K^a-1))-dm*M-u*M;                               % M popultion
  dydt(4) = rp*P*(N1+N2+N3)- sp*(N1+N2+N3)*P + u*M - dm*P;                                                      % P population
  dydt(5) = lambda*rho*N3 + alpha3*N2*N3 - alpha4*N3 + alpha5*N1*N3 - alpha6*N3 + gamma*N2 - delta3*sp*N3*P;    % N3 population

  % Alternative options for N1 and N2 conformity as in SI

  %dydt(1) = r1*N1 - alpha1*N1*N2.*(1+D.*(1-n2).*(2*n2-1)) + alpha2*N2;   % N1 population
  %dydt(2) = rho*N2 + alpha1*N1*N2.*(1+D.*(1-n2).*(2*n2-1)) -alpha2*N2;   % N2 population
  
  dydt = dydt';

end