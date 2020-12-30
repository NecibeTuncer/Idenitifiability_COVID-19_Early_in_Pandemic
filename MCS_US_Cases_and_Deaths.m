clear all
close all
clc

global tforward tmeasure_cases initial_cond tmeasure_deaths

numiter = 1000;

%true_params = [0.0176105528770787,12.1812505304075,0.00437893442353512,0.152401298752873,0.973517961499494,0.127611395979819,531.820564657953,0.00111581478655476,2290.21298128071,0.00939933012084721];
true_params = [0.0672864291950681,0.00247929868857355,0.00418910589265833,0.159033603989022,0.999725815859592,0.00912452195013469,1429.39081369128,0.114694445495184,1576.42617609149,0.000321233733707682];
%params = [ beta q k rho gamma nu alpha r];
dt = 0.1; 

tdata_cases = 1:1:53;
tdata_deaths = 2:1:53;
tdata_deaths = [0 tdata_deaths];
tforward = (0:dt:53); 
tmeasure_deaths = 21:1/dt:length(tforward);
tmeasure_deaths = [1 tmeasure_deaths];
tmeasure_cases = 11:1/dt:length(tforward);

initial_cond = [50000000  2500 1000 63 20 0];
lb = [0 0 0.0 0.0 0 0 0 0 0 0];  
ub = [1e+5 1e+5 0.35 1 1 1 1e+5 1 1e+5 1];

[~, y_r]= ode15s(@(t,y)Model_Corona(y,true_params),tforward,initial_cond);


Cases_trp = true_params(4)*true_params(3)*y_r(tmeasure_cases(:),2);
 
Deaths_trp = true_params(6)*y_r(tmeasure_deaths(:),5);


X = zeros(length(true_params),numiter);


noiselevel = [0, 0.01, 0.05, 0.1, 0.2];
total_ARE =  zeros(length(noiselevel), length(true_params));

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
 
   X(:,i) = k';
end

arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
        
    end
    
    total_ARE(noisei,:) = arescore;
    
    % ts = tinv([0.05  0.95],length(X(1,:))-1);
    % CI = mean(X(3,:)) + ts*(std(X(3,:))/sqrt(length(X(3,:))))
end

save('MCS_US_Cases_And_Deaths')
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

beta = params(1);
q = params(2);
kk = params(3);
rho = params(4);
gamma_n = params(5);
nu = params(6);
alpha = params(7);
r = params(8);
q1 = params(9);
gamma_s = params(10);

S  = y(1);
E  = y(2);
In = y(3);
Is = y(4);
J  = y(5);
R  = y(6);

N  = y(1) + y(2) + y(3) + y(4) + y(6);


dy(1) =  - beta*(q*E + q1*In +Is).*S./N ;
dy(2) =  beta*(q*E + q1*In +Is).*S./N   - kk*E;
dy(3) = (1-rho)*kk*E - gamma_n*In;
dy(4) = rho*kk*E - (alpha + gamma_s)*Is;
dy(5) = alpha*Is - (r + nu)*J;
dy(6) = gamma_n*In + r*J + gamma_s*Is;

end