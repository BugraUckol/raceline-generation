%% Prep
clc, clear%, close all
% Clear is especially important since the opti object should not be 
%   given constraints from the past

%% Casadi Imports
addpath("/Users/bugrauckol/Documents/share/casadi-3")
% addpath("C:\Program Files\casadi-3.6.7-windows64-matlab2018b")
import casadi.*

%% Import path properties
path = load('three_d_infinity.mat');

%% System Model
%{

X = [t, ey, ez, e_psi, e_the, e_phi]
    - Time
    - Error in y direction of the Frenet Frame
    - Error in z direction of the Frenet Frame
    - Yaw of the Body Frame wrt. Frenet Frame
    - Pitch of the Body Frame wrt. Frenet Frame
    - Roll of the Body Frame wrt. Frenet Frame
        ________
-[X]-->|        |
-[p]-->| System |---[X]->
-[q]-->|________|

%}

%% Constants
v0 = 1.0; % Not a state for constant velocity model
pq_lim = 1;

%% Initial conditions
t0 = 0;
e_y0 = 0.0;
e_z0 = 0.0;
e_psi = 0.0;
e_the = 0.0;
e_phi = 0.0;

%% Setting Optimization Problem
size_vec = size(s_arr);
N = size_vec(1) - 1;

opti = casadi.Opti();

X = opti.variable(6, N+1); % state trajectory in path frame
t = X(1,:);
en = X(2,:);
eb = X(3,:);
epsi = X(4,:);
ethe = X(5,:);

U = opti.variable(2,N);   %steering

% Cost Function
opti.minimize(1.0 * U(1,:) * U(1,:)' + 1.0 * U(2,:) * U(2,:)');
% Minimize angular velocity inputs p and q

% ---- dynamic constraints --------
% x' = [t, ey, ep, e_the, e_psi]
f = @(tt,een,eeb,eepsi,eethe,p,q,kappa,tau) [
    (1 - kappa * een) / (v0 * cos(eepsi) * cos(eethe));
    (1 - kappa * een) * tan(eepsi) + tau * eeb;
    (1 - kappa * een) * tan(eethe)/cos(eepsi) - tau * een;
    p - tau * sin(eepsi);
    q - kappa + tau*tan(eethe) - tau * sin(eepsi) * tan(eepsi) * tan(eethe) 
   ];