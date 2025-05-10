%% Prep
clc, clear%, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
addpath("/Users/bugrauckol/Documents/share/casadi-3")
import casadi.*
load three_d_infinity

%% Initial conditions and constants
t0 = 0;
v0 = 1.0; % Not a state for constant velocity model
e_y0 = 0.0;
e_z0 = 0.0;
e_psi = 0.0;
e_the = 0.0;

%% Vehicle Model
pq_lim = 1;

%% Setting Optimization Problem
size_vec = size(a_arr);
N = size_vec(1) - 1;

opti = casadi.Opti();

X = opti.variable(5, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
epsi = X(4,:);
ethe = X(5,:);

U = opti.variable(2,N);   %steering

% ---- objective          ---------
% opti.minimize(t(end)); % minimize time
opti.minimize(1.0 * U(1,:) * U(1,:)' + 0.0 * t(end)); % minimize steering

% ---- dynamic constraints --------
% x' = [t, ey, ep]
f = @(tt,een,eeb,eepsi,eethe,p,q,kappa,tau) [
    (1 - kappa * een) / (v0 * cos(eepsi) * cos(eethe));
    (1 - kappa * een) * tan(eepsi) + tau * eeb;
    (1 - kappa * een) * tan(eethe)/cos(eepsi) - tau * een;
    p - tau * sin(eepsi);
    q - kappa + tau*tan(eethe) - tau * sin(eepsi) * tan(eepsi) * tan(eethe) 
   ];