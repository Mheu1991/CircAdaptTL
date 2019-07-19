function PNew
% function PNew
%
% Sets backbone of structure P
% Attribute names to ArtVen, TriSeg, Chambers, Valves, Tubes, 
% Walls and Patches
% Connections between elements are defined by strings stored in
% structure field P.Net
%
% Structure is built on information in P.Net
% ArtVen represents artery-peripheral resistance-vein of organ or body part
% Chamber represents a cavity enclosed by a myocardial wall (atria)
% TriSeg represents combination of two cavities with three walls
% Bag represents passive elastic bag, encapsulating part of the
% circulation, like the pericardium
% Node: named connection point
% Cavity: volume with elastic or muscular wall, connected to a node
% Wall: muscular wall can contract, composed of 1 or more patches
% Patch: contracting part of a wall, having specific mechanical properties
% Valve: valve with inertia connects proximal to distal node, may leak
% Tube: elastic wave-guiding tube connects proximal to distal node
%
% Theo Arts, Maastricht University, April 3, 2013
% 
% Imposed changes in the calculation of vicous friction of blood flow, now based on an
% approximation that is in agreement with Womersley theory (considered gold
% standard), Maarten Heusinkveld, August 2018

global P;
P=CreateP; % sets the tree structure with all necessary fields

%======= DEFINITION OF STRUCTURE BY ELEMENTS AND CONNECTIONS ==========
Net.ArtVen.Name = {'Ca','Br','Ce','Fe','Pu'}; % carotid, brachial, celiac, femoral + Pulmonary 
Net.Chamber.Name= {'La','Ra'}; %atria
Net.TriSeg.Name = {'v'}; %ventricular unit

% Valve and Tube connections
Net.Valve.Nodes=[... % proximal and distal node of valves
    {'Vc'  ,'Ra'  };...
    {'Ra'  ,'Rv'  };...
    {'Rv'  ,'PuAr'};...
    {'PuVe','La'  };...
    {'La'  ,'Lv'  };...
    {'Lv'  ,'Ao'  }];

Net.Tube.Nodes=[... % proximal and distal node of elastic tubes
    {'Ao'  ,'BrAr'};...
    {'BrAr','CaAr'};...
    {'BrAr','CeAr'};...
    {'CeAr','FeAr'};...
    {'Vc'  ,'BrVe'};...
    {'BrVe','CaVe'};...
    {'BrVe','CeVe'};...
    {'CeVe','FeVe'}];

% Bags pressurize parts of the circulation
Net.Bag.Name   ={'Peri'     ,'Thorax'}; %pericardium, thorax
Net.Bag.Chamber={{'La','Ra'},{}      }; % enclosed chambers
Net.Bag.TriSeg ={{'v'}      ,{}      }; % enclosed TriSeg
Net.Bag.ArtVen ={{}         ,{'Pu'}  }; % enclosed ArtVen
Net.Bag.Tube   ={{}         ,{'AoBrAr','BrArCeAr','VcBrVe','BrVeCeVe'}};
Net.Bag.Bag    ={{}         ,{'Peri'}}; % pericardium inside thorax

Net.Node.Baro ={'BrAr'}; % pressure control node
Net.Wall.MultiPatch={}; % defines Wall's split in Patch's

% Transfer of structure information to P-structure
P.Net = Net;
Indexation; % mutual element relations expressed by indices

%==========================================
%============== Data filling ==============
%==========================================

% Species specific general information
P.General.q0           = 85e-6;
P.General.p0           = 12200;
P.General.tCycle       = 0.85 ;

P.General.rhob         = 1050 ;
P.General.TauAv        = 0.16; % AV depolarization delay
P.General.SaturationControl= 0  ;
P.General.AdaptFeedback    = 0.3;
rho                  = 1.05e3;  % blood density             [kg/m3]
eta                  = 3e-3;    % dynamic blood viscosity   [kg/m s]
P.General.Viscosity  = eta/rho; % blood viscosity            [m^2 s]
P.General.Profile    = 'Approximate'; %Flat % choice of velocity profile in blood vessels (choose 'Approximate' in case blood flow in medium to small (e.g. diam ~1-5mm are considered)

% Crude estimates of maximum pressures
pLv= 15000  ; % peak Lv pressure
pRv= 3800   ; % peak Rv pressure
pLa= 1400   ; % mean pLa
pRa= 560    ; % mean pRa
pPu= 1900   ; % mean pulmonary artery pressure
p0 = P.General.p0;

% General settings
P.General.tCycleRest = P.General.tCycle;
P.General.FacpControl= 1    ;
P.General.Dt         = 0.002;
P.General.ScaleVqY   = [1e-5,1e-4,1e-1];
P.General.tEnd       = 1.0;
P.General.TimeFac    = 1  ;

%======= Adaptation setpoints

AdaptationParameters
SarcomereProperties
MakeHeart(pLv,pRv,pLa,pRa)
ArtVenParameters(p0,pPu,pLa,pRa)
TubeParameters
BagParameters
ValveParameters

P.t=0;
P2SVar;

save P P

end

% ================= Auxilary functions ====================
function AdaptationParameters

global P;

% ArtVen.Adapt and Tube.Adapt
Put({'ArtVen','Adapt'},'WallStress','All',[1;1]*500e3);
Put({'ArtVen','Adapt'},'vFlowMean' ,'All',[1;1]*0.17 );
Put({'ArtVen','Adapt'},'vImpact'   ,'All',[1;1]*3.0  );
Put({'Tube'  ,'Adapt'},'WallStress','All',500e3);
Put({'Tube'  ,'Adapt'},'vFlowMean' ,'All',0.17 );
Put({'Tube'  ,'Adapt'},'vImpact'   ,'All',3.0  );

% Patch/Sarcomere
PatchA={};
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)];
    end
end
% Ventricle=default
Put({'Patch','Adapt'},'SfPasMax','All' ,6000 );
Put({'Patch','Adapt'},'SfPasAct','All' ,4800 );
Put({'Patch','Adapt'},'FacSfAct','All' ,0.61 );
Put({'Patch','Adapt'},'LsPasAct','All' ,2.1  );
Put({'Patch','Adapt'},'LsPasActT','All',2.1 ); % This line was lacking in the CircAdapt0219a version -> caused adaptation to crash
Put({'Patch','Adapt'},'SfPasMax',PatchA,30000);
Put({'Patch','Adapt'},'FacSfAct',PatchA,0.35 );

end

function SarcomereProperties
global P;
PatchA={};
for iP=1:P.Patch.n
    if P.Patch.Name{iP}(2)=='a'
        PatchA=[PatchA,P.Patch.Name(iP)];
    end
end

% Ventricular: default
Put('Patch','Lsi'             ,'All' , 1.9      );
Put('Patch','C'               ,'All' , 0.001    );
Put('Patch','Depolarization'  ,'All' , [0;0]-P.General.tCycle);
Put('Patch','LsRef'           ,'All' , 2.0      );
Put('Patch','Ls0Pas'          ,'All' , 1.8      );
Put('Patch','dLsPas'          ,'All' , 0.6      );
Put('Patch','SfPas'           ,'All' , 22000    );
Put('Patch','Lsi0Act'         ,'All' , 1.51     );
Put('Patch','LenSeriesElement','All' , 0.04     );
Put('Patch','SfAct'           ,'All' , 84000    );% lowered P114
Put('Patch','vMax'            ,'All' , 7.0      );
Put('Patch','TR'              ,'All' , 0.25     );
Put('Patch','TD'              ,'All' , 0.25     );
Put('Patch','TimeAct'         ,'All' , 0.42     );
Put('Patch','CRest'           ,'All' , 0.0      );

%Atrial: non-default
Put('Patch','SfPas'           ,PatchA, 50000    );
Put('Patch','SfAct'           ,PatchA, 59000    );% lowered P114
Put('Patch','TimeAct'         ,PatchA, 0.15     );
Put('Patch','TR'              ,PatchA, 0.4      );
Put('Patch','TD'              ,PatchA, 0.4      );

end
%=============================

function MakeHeart(pLv,pRv,pLa,pRa)
global P;
PatchA={'La1','Ra1'};
PatchV={'Lv1','Sv1','Rv1'};

% Geometry of heart walls
VStroke = P.General.q0*P.General.tCycle;
A0      = (580*VStroke^2)^(1/3);
Put('Patch','VWall',PatchA, VStroke*[28,32  ].*[pLa,pRa    ]./Get('Patch','SfAct',PatchA));
Put('Patch','VWall',PatchV, VStroke*[11,5,22].*[pLv,pLv,pRv]./Get('Patch','SfAct',PatchV));
Put('Patch','AmRef',PatchA, 0.53           *A0);
Put('Patch','AmRef',PatchV,[0.67,0.33,0.86]*A0);

% Cavity starting values of volume state variables V
% Heart
Heart={'La','Ra','Lv','Rv'};
Put('Cavity','V','All', 0); %initialization
Put('Cavity','V',Heart, VStroke*[1.06 0.73 1.69 1.52]);
%TriSeg
Put('TriSeg','V','v',0.52*VStroke      );
Put('TriSeg','Y','v',0.80*VStroke^(1/3));
P.TriSeg.Tau = 2.5*P.General.Dt; % lowpass for V and Y pre-estimate

end

function ArtVenParameters(p0,pPu,pLa,pRa)
global P
Sy = {'Ca','Br','Ce','Fe'};
Pu = 'Pu';
q0 = P.General.q0;

Put('ArtVen','q0AV' ,Sy, [0.15,0.12,0.57,0.16]*q0);
Put('ArtVen','q0AV' ,Pu, q0);
%   peripheral flow distribution over ArtVen vessel beds
Put('ArtVen','kAV'  ,Sy, 1);
Put('ArtVen','kAV'  ,Pu, 2);
Put('ArtVen','p0'   ,Sy,[p0 ; pRa]);
Put('ArtVen','p0'   ,Pu,[pPu; pLa]);
Put('ArtVen','k'    ,'All',[12;12]); %estimate
Put('ArtVen','Len'  ,'All',6.0*P.ArtVen.q0AV.^(1/3)); %estimate
Put('ArtVen','Len'  ,{'Br','Fe'},[0.3,0.5]); %arm and leg are long
P.ArtVen.A0=([1;1]*P.ArtVen.q0AV)./P.ArtVen.Adapt.vFlowMean;
P.ArtVen.AWall= diag([0.20;0.10])*P.ArtVen.A0;
P.ArtVen.p0AV = P.ArtVen.p0(1,:)-P.ArtVen.p0(2,:);

% Cavities of Artven
iCavArt= P.ArtVen.iCavity;
iCavVen= iCavArt+1;
V      = P.ArtVen.A0 .* P.ArtVen.Len([1,1],:);
P.Cavity.V(iCavArt)=V(1,:);
P.Cavity.V(iCavVen)=V(2,:);

end

function TubeParameters
global P
% Flow distribution calculated from ArtVen flows and shunt flows
FlowDistribution;  % mean flow distribution in ArtVen, Tubes, Valves
P.Tube.A0=abs(P.Tube.q./P.Tube.Adapt.vFlowMean); % Tube cross-section

% Sequence TubeNames:
Art= {'AoBrAr','BrArCaAr','BrArCeAr','CeArFeAr'};
Ven= {'VcBrVe','BrVeCaVe','BrVeCeVe','CeVeFeVe'};
All=[Art,Ven];
p0 = P.General.p0; % arterial pressure
q0 = P.General.q0;
Aux= Get('ArtVen','p0','Br');
pRa= Aux(2); % venous pressure

% pressure distribution in tubes
Put('Tube','p0',Art,p0 );
Put('Tube','p0',Ven,pRa);
% additional tube properties
Put('Tube','k'      ,All, repmat([8,11,11,15],[1,2]) );
Put('Tube','Len'    ,All, repmat([7,12,23,15]/100,[1,2]));
Put('Tube','AWall'  ,All, P.Tube.A0.* ...
    (12*P.Tube.p0./P.Tube.Adapt.WallStress+P.Tube.Adapt.vImpact*0.02));

% Volume and flow
Put('Tube','V','All', P.Tube.Len.*P.Tube.A0);
Put('Tube','q','All', P.Tube.A0.*P.Tube.Adapt.vFlowMean);

% Memory allocation for Wave delay lines, needed for Tube function
Row=zeros(1 ,P.Tube.n);
P.Tube.pL = P.Tube.p0;% Left wave, non-delayed, circular storage
P.Tube.pR = P.Tube.p0;
P.Tube.pP = P.Tube.p0;% Delayed signal
P.Tube.pD = P.Tube.p0;
P.Tube.DiL= Row+1;% delay in time step units
P.Tube.DiR= Row+1;
end

function BagParameters
global P;
Heart={'La','Ra','Lv','Rv'};
% Pericardium, Thorax Bags
P.Bag.k     = [10,10];
P.Bag.pAdapt= [200,50]; % transmural bag pressure pressure
%Estmate heart volume
VWall = sum(Get('Patch','VWall','All'));
VCav  = sum(Get('Cavity','V',Heart));
VHeart= VCav+VWall;
P.Bag.VRef= [1.4*VHeart,2.7*VHeart]; % thorax volume not really used
end

function ValveParameters

A=Get('Tube','A0','AoBrAr'); % Aortic cross-section

Put('Valve','q'    ,'All',0.0);
Put('Valve','AOpen','All',A  );
Put('Valve','ALeak','All',A*1e-6);
Put('Valve','Len'  ,'All',sqrt(A));
% Mitral and Tricuspid valve are larger
Put('Valve','AOpen',{'RaRv','LaLv'},...
    1.5* Get('Valve','AOpen',{'RaRv','LaLv'}) );
% vene-atrial orifices are inertias without valve leaflets
Put('Valve','ALeak',{'VcRa','PuVeLa'},...
    Get('Valve','AOpen',{'VcRa','PuVeLa'}) );
% L-R shunting valves are closed

% Wall: AmDead (non-contractile area)
ALv=sum(Get('Valve','AOpen',{'LvAo'  ,'LaLv'}));
ARv=sum(Get('Valve','AOpen',{'RvPuAr','RaRv'}));
ALa=sum(Get('Valve','AOpen',{'PuVeLa','LaLv'}));
ARa=sum(Get('Valve','AOpen',{'VcRa'  ,'RaRv'}));
Put('Wall','AmDead',{'Lv','Rv','La','Ra'},[ALv,ARv,ALa,ARa]);
end

function FlowDistribution
% calculates flow distribution through ArtVens, Valves and Tubes
% Flow ~=0 are assumed to be known. All other flows are unknown. If
% a sufficient number of flows is known, equations on steady state flow
% distribution are solved in a least squares sense

global P

nNode=P.Node.n;
nValve=P.Valve.n;
nTube=P.Tube.n;
nArtVen=P.ArtVen.n;

Ma =zeros(nNode,nArtVen);
Mt =zeros(nNode,nTube);
Mv =zeros(nNode,nValve);

for ia=1:nArtVen
    iP=P.ArtVen.iNode(ia);
    iD=iP+1;
    Ma(iP,ia)=Ma(iP,ia)-1;
    Ma(iD,ia)=Ma(iD,ia)+1;
end
for ia=1:nTube
    iP=P.Tube.iNodeProx(ia);
    iD=P.Tube.iNodeDist(ia);
    Mt(iP,ia)=Mt(iP,ia)-1;
    Mt(iD,ia)=Mt(iD,ia)+1;
end
for ia=1:nValve
    iP=P.Valve.iNodeProx(ia);
    iD=P.Valve.iNodeDist(ia);
    Mv(iP,ia)=Mv(iP,ia)-1;
    Mv(iD,ia)=Mv(iD,ia)+1;
end
M=[Ma,Mt,Mv];
Rga= 1:nArtVen; Rgt= nArtVen+(1:nTube); Rgv= nArtVen+nTube+(1:nValve);
nM = size(M,2);
q=[P.ArtVen.q0AV,P.Tube.q(1,:),P.Valve.q(1,:)];
Rg0=find(q~=0); %known flows q
Rg1=setdiff(1:nM,Rg0); % unknown flows q
q0=q(Rg0);
q1=-pinv(M(:,Rg1))*M(:,Rg0)*q0';
q(Rg1)=q1;
P.ArtVen.q0AV=q(Rga);
P.Tube.q=q(Rgt);
P.Valve.q=q(Rgv);
end


function P=CreateP
% Creates empty structure P with fields
P=[];

FieldsP={
    'General'
    'ArtVen'
    'Chamber'
    'TriSeg'
    'Valve'
    'Tube'
    'Node'
    'Cavity'
    'Wall'
    'Patch'
    'Bag'
    'Net'
    'SVar'
    'SVarDot'
    't'
    'tDot'
    };

FieldsGeneral={
    'q0'
    'p0'
    'tCycle'
    'tCycleRest'
    'DtSimulation'
    'tStart'
    'tEnd'
    'Dt'
    'ScaleVqY'
    'FacpControl'
    'TimeFac'
    'TauAv'
    'rhob'
    'AdaptFunction'
    'Fast'
    'In'
    'Out'
    };

FieldsArtVen={
    'Name'
    'n'
    'iCavity'
    'iNode'
    'k'
    'Len'
    'A0'
    'p0'
    'AWall'
    'p0AV'
    'q0AV'
    'kAV'
    'q'
    'qProx'
    'qDist'
    'Adapt'
    'ArtVen2NodeArt'
    'ArtVen2NodeVen'
    };

FieldsChamber={
    'Name'
    'n'
    'iCavity'
    'iWall'
    };

FieldsTriSeg={
    'Name'
    'n'
    'iCavity'
    'iWall'
    'V'
    'Y'
    'VDot'
    'YDot'
    'VS'
    'YS'
    'Tau'
    };

FieldsCavity={
    'Name'
    'n'
    'iNode'
    'V'
    'VDot'
    'A'
    'Z'
    'p'
    'pTrans'
    'Cavity2Node'
    };

FieldsWall={
    'Name'
    'n'
    'nPatch'
    'iPatch'
    'VWall'
    'Am0'
    'AmDead'
    'DADT'
    'T'
    'Cm'
    'Am'
    'pTrans'
    };

FieldsPatch={
    'Name'
    'n'
    'Lsi'
    'C'
    'LsiDot'
    'CDot'
    'Depolarization'
    'LsRef'
    'SfPas'
    'Ls0Pas'
    'dLsPas'
    'SfAct'
    'Lsi0Act'
    'LenSeriesElement'
    'vMax'
    'TimeAct'
    'TR'
    'TD'
    'CRest'
    'VWall'
    'AmRef'
    'Ef'
    'Ls'
    'SfEcm'
    'SfTit'
    'Sf'
    'DSfDEf'
    'T'
    'DADT'
    'Am0'
    'Am'
    'TauRefrac'
    'DepPath'
    'iPace'
    'Adapt'
    };

FieldsNode={
    'Name'
    'n'
    'iBaro'
    'q'
    'p'
    'Y'
    'A'
    };

FieldsBag={
    'Name'
    'n'
    'iCavity'
    'iWall'
    'iTube'
    'iBag'
    'OK'
    'VRef'
    'k'
    'pAdapt'
    'V'
    'pTrans'
    'p'
    'VC'
    'VW'
    'VT'
    'VB'
    };

FieldsValve={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'AOpen'
    'ALeak'
    'Len'
    'q'
    'qDot'
    'AvValves'
    'AvWalls'
    'Valve2NodeProx'
    'Valve2NodeDist'
    };

FieldsTube={
    'Name'
    'n'
    'iNodeProx'
    'iNodeDist'
    'k'
    'Len'
    'TauAtt'
    'A0'
    'p0'
    'AWall'
    'V'
    'VDot'
    'q'
    'A'
    'pTrans'
    'p'
    'pProx'
    'pDist'
    'qProx'
    'qDist'
    'ZR'
    'ZL'
    'pL'
    'pR'
    'pP'
    'pD'
    'DiL'
    'DiR'
    'cL'
    'cR'
    'Tube2NodeProx'
    'Tube2NodeDist'
    'Adapt'
    };

FieldsNet={
    'ArtVen'
    'Chamber'
    'TriSeg'
    'Node'
    'Valve'
    'Tube'
    'Bag'
    'Wall'
    'Depolarization'
    };

Empty10=ones(1,0);
for i=1:length(FieldsP)
    fld=FieldsP{i};
    P.(fld)=[];
    FldStr=['Fields',fld];
    if exist(FldStr)
        SubFldStr=eval(FldStr);
        for j=1:length(SubFldStr)
            P.(fld).(SubFldStr{j})=Empty10;
        end
    end
end

end

