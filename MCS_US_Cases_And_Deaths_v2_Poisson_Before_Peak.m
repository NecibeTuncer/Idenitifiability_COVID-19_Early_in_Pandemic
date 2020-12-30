clear all
close all
clc

global tforward tmeasure_cases initial_cond tmeasure_deaths

numiter = 1000;

% true_params = [0.0718083403572489,0.00459800674201816,0.00418211805318638,0.130568075558499,0.999998302510065,0.00454681393625922,266.192220147055,1.18387211195392e-06,1562.08224106316,0.231967881518898]
% true_params(2) = true_params(1)*true_params(2);
% true_params(9) = true_params(1)*true_params(9);

true_params = [0.0718,0.0003,0.0042,0.1306,1,0.0045,266.1922,1.1e-6,112.1705,0.232];
%params = [ beta q k rho gamma nu alpha r];
dt = 0.1; 

tdata_cases = 1:1:28;
tdata_deaths = 2:1:28;
tdata_deaths = [0 tdata_deaths];
tforward = (0:dt:28); 
tmeasure_deaths = 21:1/dt:length(tforward);
tmeasure_deaths = [1 tmeasure_deaths];
tmeasure_cases = 11:1/dt:length(tforward);

initial_cond = [50000000  2500 1000 63 20 0];
lb = [0 0 0.0 0.0 0 0 0 0 0 0];  
%ub = [1e+5 1e+5 0.35 1 1 1 1e+5 1 1e+5 1];
 S0N0ratio  = initial_cond(1)/(initial_cond(1)+initial_cond(2)+initial_cond(3)+initial_cond(4)+initial_cond(5)+initial_cond(6))

 Reproduction_Number = S0N0ratio *(true_params(2)/true_params(3) + (1-true_params(4))*true_params(9)/true_params(5) + true_params(4)*true_params(1)/(true_params(10) + true_params(7)));

[~, y_r]= ode15s(@(t,y)Model_Corona(y,true_params),tforward,initial_cond);


Cases_trp = true_params(4)*true_params(3)*y_r(tmeasure_cases(:),2);
 
Deaths_trp = true_params(6)*y_r(tmeasure_deaths(:),5);




total_ARE =  zeros(1, length(true_params));
X = zeros(length(true_params),numiter);



for i = 1:numiter
    i
    
  CoronaCases = poissrnd(Cases_trp);

  CoronaDeaths = poissrnd(Deaths_trp); 
  
   k = true_params;
   lb = zeros(size(true_params));
   
   [k,~] =  fminsearchbnd(@(k)err_in_data(k,CoronaCases,CoronaDeaths),k,lb,[]);
 
   X(:,i) = k';
end

arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
        
    end
    
    total_ARE(1,:) = round(arescore,1);
     % ts = tinv([0.05  0.95],length(X(1,:))-1);
    % CI = mean(X(3,:)) + ts*(std(X(3,:))/sqrt(length(X(3,:))))


save('MCS_US_Cases_And_Deathsv2_Poisson_Before_the_Peak')
function error_in_data = err_in_data(k,CoronaCases,CoronaDeaths) 

global tforward initial_cond tmeasure_cases tmeasure_deaths
      
       
   
 [~,y] = ode15s(@(t,y)Model_Corona(y,k),tforward,initial_cond);
 
 Cases = k(4)*k(3)*y(tmeasure_cases(:),2);
  
 Deaths = k(6)*y(tmeasure_deaths(:),5);

 w_cases = (1/((mean(CoronaCases)^2)*length(CoronaCases)))*10;
 w_deaths = 1/((mean(CoronaDeaths)^2)*length(CoronaDeaths));
 
 error_in_data = sum((Cases - CoronaCases).^2)*w_cases+...
                 sum((Deaths - CoronaDeaths).^2)*w_deaths;       
  
end

function dy = Model_Corona(y,params)


dy = zeros(6,1);

beta_s = params(1);
betaE = params(2);
kk = params(3);
rho = params(4);
gamma_n = params(5);
nu = params(6);
alpha = params(7);
r = params(8);
beta_n = params(9);
gamma_s = params(10);

S  = y(1);
E  = y(2);
In = y(3);
Is = y(4);
J  = y(5);
R  = y(6);

N  = y(1) + y(2) + y(3) + y(4) + y(6);



dy(1) =  - (betaE*E + beta_n*In + beta_s*Is).*S./N ;
dy(2) =  (betaE*E + beta_n*In + beta_s*Is).*S./N   - kk*E;
dy(3) = (1-rho)*kk*E - gamma_n*In;
dy(4) = rho*kk*E - (alpha + gamma_s)*Is;
dy(5) = alpha*Is - (r + nu)*J;
dy(6) = gamma_n*In + r*J + gamma_s*Is;

end