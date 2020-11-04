Class SIV_Script_multi.m

clear all; close all; clc;
 
p=struct;
p.N=306800000;
p.a=0.0033;
p.G=1/3;
p.A=1/14;
p.vaccine_release_date=289;
report_rate = 0.316;
hospitalization_morbidity_ratio = 274304/60800000*report_rate;
mortality_morbidity_ratio = 12469/60800000*report_rate;
 
 
Tspan = [55 450]; % beginning and ending time
 
% v(1)=S1, v(2)=S2 v(3)=I v(4)=V v(5)=M v(6)=R%
% 
Initial_Condition = [(306800000-1)*.11;(306800000-1)*.89;1;0;0;0;0]; % initial condition
 
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
vaccinated_289 = max(vsol(:,7))
vaccination_cost_289 = vaccinated_289 * 26
 
figure(7)
plot(t,vsol(:,7),'-')
 
hold on
 
p.vaccine_release_date=259;
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
vaccinated_259 = max(vsol(:,7))
vaccination_cost_259 = vaccinated_259 * 26
plot(t,vsol(:,7),'-')
 
hold on
 
p.vaccine_release_date=100;
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
vaccinated_100 = max(vsol(:,7))
vaccination_cost_100 = vaccinated_100 * 26
plot(t,vsol(:,7),'-')
 
hold on
 
p.vaccine_release_date=319;
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
vaccinated_319 = max(vsol(:,7))
vaccination_cost_319 = vaccinated_319 * 26
plot(t,vsol(:,7),'-')
 
 
p.vaccine_release_date=289;
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
figure(1)
plot(t,vsol(:,1),'-')
 
figure(2)
plot(t,vsol(:,2),'-')
 
figure(3)
plot(t,vsol(:,3),'-')
 
figure(4)
plot(t,vsol(:,4),'-')
 
figure(5)
plot(t,vsol(:,5),'-')
 
figure(6)
infected_group_289 = max(vsol(:,6))
hospitalization_289 = infected_group_289 * hospitalization_morbidity_ratio
hospital_cost_289 = hospitalization_289 * 1986 * 3
death_289 = infected_group_289 * mortality_morbidity_ratio
wage_cost_289 = infected_group_289 * 3/365*59039/2
plot(t,vsol(:,6),'-')
 
hold on
 
p.vaccine_release_date=259;
 
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
infected_group_259 = max(vsol(:,6))
hospitalization_259 = infected_group_259 * hospitalization_morbidity_ratio
hospital_cost_259 = hospitalization_259 * 1986 * 3
death_259 = infected_group_259 * mortality_morbidity_ratio
wage_cost_259 = infected_group_259*3/365*59039/2
plot(t,vsol(:,6),'-')
 
hold on
 
p.vaccine_release_date=100;
 
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
infected_group_100 =  max(vsol(:,6))
hospitalization_100 = infected_group_100 * hospitalization_morbidity_ratio
hospital_cost_100 = hospitalization_100 * 1986 * 3
death_100 = infected_group_100 * mortality_morbidity_ratio
wage_cost_100 = infected_group_100*3/365*59039/2
 
plot(t,vsol(:,6),'-')
 
hold on
p.vaccine_release_date=319;
 
[t,vsol] = ode45(@(t,v) susceptible_infected_vaccinated_multi(t,v,p),Tspan,Initial_Condition);
infected_group_319 = max(vsol(:,6))
hospitalization_319 = infected_group_319 * hospitalization_morbidity_ratio
hospital_cost_319 = hospitalization_319 * 1986 * 3
death_319 = infected_group_319 * mortality_morbidity_ratio
wage_cost_319 = infected_group_319*3/365*59039/2
plot(t,vsol(:,6),'-')
 
 
total_cost_100 = wage_cost_100 + hospital_cost_100 + vaccination_cost_100
total_cost_259 = wage_cost_259 + hospital_cost_259 + vaccination_cost_259
total_cost_289 = wage_cost_289 + hospital_cost_289 + vaccination_cost_289
total_cost_319 = wage_cost_319 + hospital_cost_319 + vaccination_cost_319



Class susceptible_infected_vaccinated_multi.m
function dvdt= susceptible_infected_vaccinated(t,v,p)
 
if t<p.vaccine_release_date
    a = 0;
else 
    a = p.a;
    
end
 
    N= p.N;
    G= p.G;
    A=p.A; 
 
    %v(1)=S1, v(2)=S2 v(3)=I v(4)=V v(5)=M v(6)=R%
    %v(7)=v2;
    
    
    dvdt = [-B(t)*(v(1)*v(3)/(N));
        -B(t)*(v(2)*v(3)/N)-a*v(2);
        (B(t)*v(2)*v(3)/N)+B(t)*(v(1)*v(3)/(N))+B(t)*v(4)*v(3)/N-G*v(3);
        a*v(2)-(A)*v(4)-(B(t)*v(4)*v(3)/N);
        A*v(4);
        G*v(3);
        a*v(2)];
    
    function Beta=B(t)
        
        Beta=.52*(1+.35*cos(2*pi*t/365));
    end
 
end

