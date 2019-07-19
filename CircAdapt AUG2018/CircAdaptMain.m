% CircAdaptMain
% Theo Arts, Tammo Delhaas, Joost Lumens
% email: t.arts@bf.unimaas.nl
%
% CircAdapt is a lumped model of the dynamic circulation. Size of heart
% and bloodvessels adapt to load, as generated by the model itself.
% All variables and parameter are stored in the structure P with file
% name 'P'. Hemodynamic data are displayed graphically
% Theo Arts, Maastricht University, Oct 13, 2012

clear
global P;

addpath(pwd);% add current directory to path

if exist('P.mat','file'); % if P.mat file is available
    FileName='P';
    PathName=[pwd,'\'];
    ShowMenu=1;
    while ShowMenu;
        c2=input('[N]ew, [R]eference, [L]oad, [C]ontinue <Enter>: ','s');
        a=[c2,'c']; c=a(1); % convert <Enter> to 'c'
        switch lower(c)
            case('l'); % load file with P-structure
                [FileName,PathName] = uigetfile('*.mat','Choose file: ');
                load([PathName,FileName]);
                c='c';ShowMenu=false;
            case('r') % load reference file with P-structure
                load('PRef');
                c='c';ShowMenu=false;
            case('n'); % new P-structure from scratch
                PNew;
                load P;
                c='c';ShowMenu=false;
            case('c'); % continuation of simulation
                load([PathName,FileName]);
                c='c';ShowMenu=false;
            otherwise
                ShowMenu=true;
        end
    end
else
    if exist('PRef.mat','file');
        c2=input('[N]ew, [R]eference, [L]oad,: ','s');
        a=[c2,'c']; c=a(1); % convert <Enter> to 'c'
        switch lower(c)
            case('l'); % load file with P-structure
                [FileName,PathName] = uigetfile('*.mat','Choose file: ');
                load([PathName,FileName]);
                c='c';ShowMenu=false;
            case('r')
                load('PRef');
                c='c';
            otherwise
                PNew;
                c='c';
        end
    else
        PNew; % Parameter initialization, some remodeling rules inclusive
        % Generates parameter structure Par and initial conditions of the variables P
        c='c';
    end
end


% Default initialization
G=P.General;
G.DtSimulation=1.5*G.tCycle; % standard duration of simulation

G.AdaptFunction='Adapt0P'; % No adaptation
G.Fast= 0; % regular beat to beat sequence

%XXXX Menu for changing hemodynamic variables and adaptation condition
OK=1; NY='NY';
while OK;
    disp(' ');
    disp(['[P]ressure                  (kPa): ',num2str(G.p0/1e3)]);
    disp(['[F]low                     (ml/s): ',num2str(G.q0*1e6)]);
    disp(['[T]ime of beat               (ms): ',num2str(G.tCycle*1e3)]);
    disp(['[D]uration simulation         (s): ',num2str(G.DtSimulation)]);
    disp(['Adapt n[O]ne [R]est,[E]xercise   : ',G.AdaptFunction]);
%     disp(['Faster steady state [Y]/[N]      : ',NY(G.Fast+1)]);% does not
% work yet++++++++++++++
    disp( '<Enter> = Continue');
    c1=input('Choose Letter <Enter>: ','s');
    switch lower(c1);
        case 'p'
            G.p0=input('Mean Arterial Pressure (kPa): ')*1e3;
        case 'f'
            G.q0=input('Systemic Flow (ml/s): ')/1e6;
        case 't'
            G.tCycle=input('Cycle Time (ms): ')*1e-3;
        case 'd'
            G.DtSimulation=input('Duration of simulation (s): ');
        case 'o'
            G.AdaptFunction='Adapt0P';
            G.DtSimulation=1.5*G.tCycle;
        case 'r'
            G.AdaptFunction='AdaptRest';
            G.DtSimulation=50*G.tCycle;
        case 'e'
            G.AdaptFunction='AdaptExc';
            G.DtSimulation=100*G.tCycle;
        case 'y'
            G.Fast=1;
            G.DtSimulation=max(30*G.tCycle,G.DtSimulation);
        case 'n'
            G.Fast=0;
        otherwise
            OK=0;
    end
end

P.General=G;
CircAdapt; %generate solution

%=== Saving State Variables and Model Parameters
save P P; %save solution as structure P in file 'P'
disp('Simulation finished successfully');

CircDisplay(0); % graphical display
% Structure P has been extended with solution for all variables as a
% function of time

