clear all
close all
clc
format long
%addpath('../Documents/Research/UtilsOptimization')

global CoronaCases CoronaDeaths tforward tmeasure_cases...
      initial_cond tmeasure_deaths


CoronaCases = [63;98;116;106;163;290;307;329;553;587;843;983;1748;2853;...
4582;5588;4825;9400;10189;11075;13355;17224;18691;19452;19913;20297;24742;...
26473;29874;32284;34196;25316;31477;33748;32210;33806;33943;30186;27608;...
26820;27088;30342;29567;32165;29057;25844;28123;26084;30156;31889;38958;...
35419;26509];
          
CoronaDeaths = [2;3;4;3;4;4;8;3;7;9;12;18;23;40;56;49;46;113;141;225;247;...
268;400;525;363;558;912;1049;974;1045;1330;1165;1496;2219;2156;2101;2226;...
2013;1715;1714;2553;2618;2176;2528;1867;1561;1939;2683;2358;2340;1957;2065;...
1157];
                    
dt = 0.1;
tdata_cases = 1:1:53;
tdata_deaths = 2:1:53;
tdata_deaths = [0 tdata_deaths];
tforward = (0:dt:53); 
tmeasure_deaths = 21:1/dt:length(tforward);
tmeasure_deaths = [1 tmeasure_deaths];
tmeasure_cases = 11:1/dt:length(tforward);
tforward_projection = 0:dt:365;


%params = [ beta q k rho gamma nu alpha r q1];
    
params = [0.0672864366335624,0.00247832863675690,0.00418910779471505,...
0.159032299140594,0.999720725613706,0.00912452737420794,1429.58768476442,...
0.114690775217238,1576.45442372945,0.000323654778013549];

params = [0.0674240060817094,0.00253400570689723,0.00417638062306131,...
0.159488452496107,1.00455283638880,0.00922730360943078,1445.21079967128,...
0.116132605909701,1582.52782427876,0.000327468968652404];

params = [0.0673840328091480,0.00266977230878318,0.00417629849415970,0.159494386616862,0.999999989344867,0.00922684835248061,1444.63625065826,0.116126485046159,1582.92374065607,0.000338858193095259];
lb = [0 0 0.0 0.0 0 0 0 0 0 0];  
ub = [1e+5 1e+5 0.35 1 1 1 1e+5 1 1e+5 1];
% lb = [0 0 0.02 0 0 0 0 0];
% ub =[1e+5 1e+5 0.16 1 1e+5 1e+5 1e+5 1e+5];
initial_cond = [50000000  2500 1000 63 20 0];

%[params,fval] =  fminsearchbnd(@err_in_data,params,lb,ub,optimset('Display','iter'))

[~, y_r] = ode15s(@(t,y)Model_Corona(y,params),tforward,initial_cond);
[~, y_p] = ode15s(@(t,y)Model_Corona(y,params),tforward_projection,initial_cond);

 
 Cases_p = params(4)*params(3)*y_r(:,2);


 Cases_end = params(4)*params(3)*y_r(tmeasure_cases(:),2);
 
 Deaths_end = params(6)*y_r(tmeasure_deaths(:),5);

 w_cases = 1/((mean(CoronaCases)^2)*length(CoronaCases));
 w_deaths = 1/((mean(CoronaDeaths)^2)*length(CoronaDeaths));
 
 error_in_Cases = sum((Cases_end - CoronaCases).^2)*w_cases;
 error_in_Deaths = sum((Deaths_end - CoronaDeaths).^2)*w_deaths; 
 
 S0N0ratio  = initial_cond(1)/(initial_cond(1)+initial_cond(2)+initial_cond(3)+initial_cond(4)+initial_cond(5)+initial_cond(6))

 Reproduction_Number = S0N0ratio *(params(1)*params(2)/params(3) + (1-params(4))*params(1)*params(9)/params(5) + params(4)*params(1)/(params(10) + params(7)))
figure(1)

plot(tforward,Cases_p,'LineWidth',2.5)
hold on 
plot(tdata_cases, CoronaCases, 'r.', 'MarkerSize',20)
H=gca;
H.LineWidth=2;
title('Corona Cases','fontweight','normal','fontsize',18)
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 10 30 53])
xticklabels({'March 3','March 13',' April 3',' April 26'})
yticks([1000 5000 10000 20000 30000 40000])
yticklabels({'1000', '5000', '10000', '20000', '30000', '40000'})
% figure(2)
% 
% plot(tforward_projection,params(4)*params(3)*y_p(:,2),'LineWidth',2.5)
% hold on 
% plot(tdata_cases, CoronaCases, 'r.', 'MarkerSize',20)
% title('Corona Cases')



figure(2)
plot(tforward,params(6)*y_r(:,5),'LineWidth',2.5)
hold on 
plot(tdata_deaths, CoronaDeaths, 'r.', 'MarkerSize',20)
H=gca;
H.LineWidth=2;
title('Corona Deaths','fontweight','normal','fontsize',18)
xlabel('Time (days)','fontweight','normal','fontsize',18)
xticks([0 10 30 53])
xticklabels({'March 3','March 13',' April 3',' April 26'})
% 
% figure(4)
% plot(tforward_projection,params(10)*y_p(:,5),'LineWidth',2.5)
% hold on 
% plot(tdata_deaths, CoronaDeaths, 'r.', 'MarkerSize',20)
% title('Corona Deaths')
% 
% figure(5)
% plot(tforward_projection,initial_cond(1) - y_p(:,1),'LineWidth',2.5)
% title('cumulative incidences')
% figure(6)
% plot(tforward_projection, y_p(:,1),'LineWidth',2.5)
% 
% figure(7)
% plot(tforward_projection, y_p(:,7),'LineWidth',2.5)
% title('cumulative counted incidences')

function error_in_data = err_in_data(k) %lsqcurvefit

global tforward initial_cond tmeasure_cases ...
        CoronaCases  CoronaDeaths tmeasure_deaths
   
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