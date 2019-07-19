% Adaptation Framework 
% Exc-Rest-Exc-Rest-Exc-Rest
% M Heusinkveld

clear all;
close all;
clc;

% Initial P-struct (non-adapted)
load('PRef.mat')

% Perform Exercise adaptation (Exc1)
SimName                 = 'Exc1';
P.General.tCycle        = 0.426;
P.General.q0            = 255e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptExc';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
save P P
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Rest adaptation (Rest1)
load([SimName,'.mat']);
SimName    = 'Rest1';
P.General.tCycle        = 0.825;
P.General.q0            = 85e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptRest';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Exercise adaptation (Exc2)
load([SimName,'.mat']);
SimName                 = 'Exc2';
P.General.tCycle        = 0.426;
P.General.q0            = 255e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptExc';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Rest adaptation (Rest2)
load([SimName,'.mat']);
SimName    = 'Rest2';
P.General.tCycle        = 0.825;
P.General.q0            = 85e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptRest';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Exercise adaptation (Exc3)
load([SimName,'.mat']);
SimName                 = 'Exc3';
P.General.tCycle        = 0.426;
P.General.q0            = 255e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptExc';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Rest adaptation (Rest3)
load([SimName,'.mat']);
SimName    = 'Rest3';
P.General.tCycle        = 0.825;
P.General.q0            = 85e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptRest';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Exercise adaptation (Exc4)
load([SimName,'.mat']);
SimName                 = 'Exc4';
P.General.tCycle        = 0.426;
P.General.q0            = 255e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptExc';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');

% Perform Rest adaptation (Rest3)
load([SimName,'.mat']);
SimName    = 'Rest4';
P.General.tCycle        = 0.825;
P.General.q0            = 255e-6; %(3*q0_baseline)
P.General.AdaptFunction = 'AdaptRest';
P.General.DtSimulation  = 150*P.General.tCycle; %(convergence is slow)
CircAdapt;
CircDisplay
P.General.DtSimulation  = 2.1 *P.General.tCycle;
save P P
CircAdapt;
CircDisplay;
save(SimName,'P');



