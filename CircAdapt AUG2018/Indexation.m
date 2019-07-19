function Indexation
% function Indexation
% Sets backbone of structure P
% Attribute names to ArtVen, TriSeg, Chambers, Valves, Tubes, 
% Walls and Patches
% Connections between elements are defined by strings stored in
% structure field P.Net
%
% The P-Structure is built on information in P.Net
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
% Tube: elastic tube connects proximal to distal node
%
% Theo Arts, Maastricht University, Feb 19, 2014

global P;

Net=P.Net; % contains complete information about element composition 
% and mutual connections of Heart and Circulation

%==== All about Name-strings of elements =============

% ArtVen->Cavity
ArtVenName   = Net.ArtVen.Name; % getting ArtVen names
ArtVen2Cavity= strcat(ArtVenName,{'Ar'}); % arterial coupled cavity names
Aux=[strcat(ArtVenName,{'Ar'});strcat(ArtVenName,{'Ve'})];
CavityName= Aux(:)'; % names of arterial and venous coupled cavities

% Chamber->[Cavity,Wall]
ChamberName    = Net.Chamber.Name; % getting chamber names
Chamber2Cavity = ChamberName; % Naming of related cavities
CavityName     = [CavityName,Chamber2Cavity]; % append chamber related cavities
Chamber2Wall   = Chamber2Cavity; % chamber related cavity names
WallName       = Chamber2Wall; % chamber related walls

% TriSeg->[Cavity,Wall]
TriSegName   = Net.TriSeg.Name; % getting triseg names
TriSeg2Cavity= strcat({'L'},TriSegName); % triseg related cavities
Aux= [strcat({'L'},TriSegName);strcat({'R'},TriSegName)];
CavityName   = [CavityName,Aux(:)']; % Append L and R cavities to cavities
TriSeg2Wall  = strcat({'L'},TriSegName);
Aux=[strcat({'L'},TriSegName);...
    strcat({'S'},TriSegName);...
    strcat({'R'},TriSegName)];
WallName     = [WallName,Aux(:)']; % Name of L, S and R wall

% Cavity->Node
Cavity2Node  = CavityName; % cavity related nodes with the same name
NodeName     = Cavity2Node;

% Valve->NodeProx/Dist
ValveNode     = Net.Valve.Nodes;% getting node names prox and dist of valve
Valve2NodeProx= ValveNode(:,1)';
Valve2NodeDist= ValveNode(:,2)';
nr=size(ValveNode,1);
ValveName=cell(1,nr);
for i=1:nr % name valves to concatenation prox and dist node
    ValveName(i)=strcat(Valve2NodeProx(i),Valve2NodeDist(i));
end
dNodeName= setdiff([Valve2NodeProx,Valve2NodeDist],NodeName);
% remove double names 
NodeName = [NodeName,dNodeName];

% Tube->NodeProx/Dist
TubeNode     = Net.Tube.Nodes;% getting node names prox and dist of tubes
Tube2NodeProx= TubeNode(:,1)';
Tube2NodeDist= TubeNode(:,2)';
nr=size(TubeNode,1);
TubeName=cell(1,nr);
for i=1:nr % name tubes to concatenation prox and dist node
    TubeName(i)=strcat(Tube2NodeProx(i),Tube2NodeDist(i));
end
dNodeName= setdiff([Tube2NodeProx,Tube2NodeDist],NodeName);
% remove double names 
NodeName = [NodeName,dNodeName];

% Wall->Patch, Wall.nPatch
Wall2Patch= strcat(WallName,'1'); % Patches are named to wall and numbered
PatchName = Wall2Patch;
WallnPatch=ones(1,length(WallName)); % 1st patch of each wall
% Walls composed with MultiPatch
MultiPatchWall= Net.Wall.MultiPatch; % search multi-patched walls
nMultiPatch   = length(MultiPatchWall);
for i=1:nMultiPatch
    NameMulti={};
    Str=MultiPatchWall{i};
    LtrStr   = isletter(Str);%logical letter positions
    NameStr  = Str(LtrStr);
    NumberStr= str2double(Str(~LtrStr));
    iWall    = strcmp( NameStr     ,WallName )*(1:numel(WallName))';
    iPatch   = strcmp([NameStr,'1'],PatchName)*(1:numel(PatchName))';
    for k=2:NumberStr
        NameMulti=[NameMulti,[NameStr,num2str(k)]];
    end
    PatchName=[PatchName(1:iPatch),NameMulti,PatchName(iPatch+1:end)];
    WallnPatch(iWall)= NumberStr; % store number of patches in wall record
end

BagName    = Net.Bag.Name;

%==== END All about Name strings =============

%==== indexations ==========
% Naming and counting of elements, determined by name strings
P.ArtVen = Naming('ArtVen' ,ArtVenName );
P.Chamber= Naming('Chamber',ChamberName);
P.TriSeg = Naming('TriSeg' ,TriSegName );
P.Valve  = Naming('Valve'  ,ValveName  );
P.Tube   = Naming('Tube'   ,TubeName   );
P.Node   = Naming('Node'   ,NodeName   );
P.Cavity = Naming('Cavity' ,CavityName );
P.Wall   = Naming('Wall'   ,WallName   );
P.Patch  = Naming('Patch'  ,PatchName  );

% indices determine mutual relations between elements
P.ArtVen.iCavity = Get('Cavity','Index',ArtVen2Cavity );
P.Chamber.iCavity= Get('Cavity','Index',Chamber2Cavity);
P.Chamber.iWall  = Get('Wall'  ,'Index',Chamber2Wall  );
P.TriSeg.iCavity = Get('Cavity','Index',TriSeg2Cavity );
P.TriSeg.iWall   = Get('Wall'  ,'Index',TriSeg2Wall   );
P.Valve.iNodeProx= Get('Node'  ,'Index',Valve2NodeProx);
P.Valve.iNodeDist= Get('Node'  ,'Index',Valve2NodeDist);
P.Tube.iNodeProx = Get('Node'  ,'Index',Tube2NodeProx );
P.Tube.iNodeDist = Get('Node'  ,'Index',Tube2NodeDist );
P.Cavity.iNode   = Get('Node'  ,'Index',Cavity2Node   );
P.Wall.nPatch    = WallnPatch;
P.Wall.iPatch    = Get('Patch' ,'Index',Wall2Patch    );

% additional through-indexation (for convenience)
P.ArtVen.iNode   = P.Cavity.iNode(P.ArtVen.iCavity);

% identify element content of bags by indices
P.Bag.Name   = BagName;
nBag         = length(P.Bag.Name);
P.Bag.n      = nBag;
P.Bag.OK     = zeros(1,P.Bag.n); % keeps track of nesting of Bags
P.Bag.iCavity= [];
P.Bag.iWall  = [];
P.Bag.iTube  = [];
P.Bag.iBag   = [];
for i=1:P.Bag.n
    BagIndex(i); % determines indices of bag-enclosed parts
end

% Baroreceptor location
P.Node.iBaro = Get('Node','Index',Net.Node.Baro);
% P.Patch.iPace= Get('Patch','Index',Net.Depolarization.Pace);
end


% ======== AUXILARY FUNCTIONS ============

function Element= Naming(ElementType,ElementName)
% Naming
global P;
Element=P.(ElementType);
Element.Name= ElementName;
Element.n   = length(Element.Name);
end

function BagIndex(iBag)
% function BagIndex(iBag)
% determines from indices of Chambers, TriSegs, ArtVens -> the indices of
% related Tubes, Walls and Cavities, enclosed by each Bag.

global P;
if P.Bag.OK(iBag)==1 % if indices already determined, skip this
    return; 
end;
iCh  = Get('Chamber','Index',P.Net.Bag.Chamber{iBag});
iTr  = Get('TriSeg' ,'Index',P.Net.Bag.TriSeg{iBag} );
iAv  = Get('ArtVen' ,'Index',P.Net.Bag.ArtVen{iBag} );
iCCh = P.Chamber.iCavity(iCh);
iCTr1= P.TriSeg.iCavity(iTr) ; iCTr=[iCTr1,iCTr1+1];
iCAv1= P.ArtVen.iCavity(iAv) ; iCAv=[iCAv1,iCAv1+1];
P.Bag.iCavity{iBag}= [iCCh,iCTr,iCAv];

iWCh = P.Chamber.iWall(iCh)  ;
iWTr1= P.TriSeg.iWall(iTr)   ;
iWTr =[iWTr1,iWTr1+1,iWTr1+2];
P.Bag.iWall{iBag}=[iWCh,iWTr];

P.Bag.iTube{iBag}= NonZero(Get('Tube','Index',P.Net.Bag.Tube{iBag}));
jBag             = NonZero(Get('Bag' ,'Index',P.Net.Bag.Bag{iBag} ));
P.Bag.iBag{iBag} = jBag;
for j=1:length(jBag)
    jB=jBag(j);
    BagIndex(jB);
    P.Bag.iCavity{iBag}= [P.Bag.iCavity{iBag}, P.Bag.iCavity{jB}];
    P.Bag.iWall{iBag}  = [P.Bag.iWall{iBag}  , P.Bag.iWall{jB}  ];
    P.Bag.iTube{iBag}  = [P.Bag.iTube{iBag}  , P.Bag.iTube{jB}  ];
end
P.Bag.OK(iBag)=1;
end

function a=NonZero(a)
a=a(logical(a~=0));
end

