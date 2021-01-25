%Hongyi Su, B.S., B.E., MSc. Student in Geophysics at LMU Munich and TUM, 30.12.2020

%This scrpt is used for two-layered(core-mentle)earth model
%Assumption:
%1. The earth is a perfect sphere(ball) 
%2. Core-mantle model, where core is the inner part and mantle is the outer part
%3. rho_mentle and rho_core are constant within its own layer, where
%rho_mentle and rho_core are density of mentle and core respectively

%Define Global Parameters 
R = 6371; % earth radius in kilometers ????????
M = 5.974*10^24; % earth mass in kg ????????
J = 8.021*10^31; % The moment of inertia of earth in kg*km^2 

%Input, Rc as propertion to R (eg: 1/3*R)
Rc = input('Rc? eg: "1/3*R" :','s');
Rc = eval(Rc); % convert format to num 
rho_mean = M/(4/3*pi*R^3) %earth mean density in kg/km^3 
Vm = 4/3*(pi)*(R^3 - Rc^3); %volume of the mantle in km^3 
Vc = 4/3*(pi)*Rc^3; %volume of the core in km^3 

%Rewirte conservation Laws(rho and J) as matrix multiplication

a = 1 - (Rc/R)^3;
b = (Rc/R)^3;
c = 2/5*Vm*(R^5-Rc^5)/(R^3-Rc^3);
d = 2/5*Vc*Rc^2;
A = [a b;c d];
want = 1/det(A)*[d -b;-c a]*[rho_mean J]'; %solve the matrix

%Example: Rc=1/3*R 
rho_mentle = want(1,1) 
rho_core = want(2,1) 
disp('unit: kg*km^2')