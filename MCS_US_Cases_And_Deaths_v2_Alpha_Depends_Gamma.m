clear all
close all
clc

global tforward tmeasure_cases initial_cond tmeasure_deaths

numiter = 1000;

% true_params = [0.00526800172928610,11.3204244573950,0.00442606656500631,0.150886999932046,0.825962513131170,0.766792178923430,0.822034532986701,13697.3246989472,0.122697074516841];
% true_params(2) = true_params(1)*true_params(2);
% true_params(8) = true_params(1)*true_params(8);
% 
% true_params = [0.00526800172928610,0.0596360156178095,0.00442606656500631,0.150886999932046,0.825962513131170,0.766792178923430,0.822034532986701,72.1575302006471,0.122697074516841];
% true_params = round(true_params,4);

true_params = [0.0053,0.0596,0.0044,0.1509,0.8260,0.7668,0.8220,72.1575,0.1227];



dt = 0.1; 

tdata_cases = 1:1:53;
tdata_deaths = 2:1:53;
tdata_deaths = [0 tdata_deaths];
tforward = (0:dt:53); 
tmeasure_deaths = 21:1/dt:length(tforward);
tmeasure_deaths = [1 tmeasure_deaths];
tmeasure_cases = 11:1/dt:length(tforward);

initial_cond = [50000000  2500 1000 63 20 0];
lb = [0 0 0.0 0.0 0 0 0 0 0 ];  
ub = [1e+5 1e+5 0.35 1 1 1 1 1e+5 1];

[~, y_r]= ode15s(@(t,y)Model_Corona(y,true_params),tforward,initial_cond);


Cases_trp = true_params(4)*true_params(3)*y_r(tmeasure_cases(:),2);
 
Deaths_trp = true_params(6)*y_r(tmeasure_deaths(:),5);





noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE =  zeros(length(noiselevel), length(true_params));
X = zeros(length(noiselevel),length(true_params),numiter);

for noisei = 1:5
    
rng default
noiselev = noiselevel(noisei)

for i = 1:numiter
    i
    
  CoronaCases = (noiselev*(Cases_trp).*randn(length(tmeasure_cases),1)) + Cases_trp;

  CoronaDeaths = (noiselev*(Deaths_trp).*randn(length(tmeasure_deaths),1)) + Deaths_trp;
 
   k = true_params;
   lb = zeros(size(true_params));
   
   [k,~] =  fminsearchbnd(@(k)err_in_data(k,CoronaCases,CoronaDeaths),k,lb,ub);
 
   X(noisei,:,i) = k';
end

arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        
        arescore(i) = 100*sum(abs(true_params(i) - X(noisei,i,:))/abs(true_params(i)))/numiter;
        
    end
    
    total_ARE(noisei,:) = round(arescore,1);
    
    % ts = tinv([0.05  0.95],length(X(1,:))-1);
    % CI = mean(X(3,:)) + ts*(std(X(3,:))/sqrt(length(X(3,:))))
end

save('MCS_US_Cases_And_Deathsv2_Alpha_Depends_Gamma')
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
r = params(7);
beta_n = params(8);
gamma_s = params(9);

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
dy(4) = rho*kk*E - (0.1765*gamma_s + gamma_s)*Is;
dy(5) = 0.1765*gamma_s*Is - (r + nu)*J;
dy(6) = gamma_n*In + r*J + gamma_s*Is;

end