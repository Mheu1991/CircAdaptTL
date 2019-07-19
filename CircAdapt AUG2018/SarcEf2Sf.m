function SarcEf2Sf
% function SarcEf2Sf
% Calculates myofiber stress Sf from myofiber strain Ef for all Patch
% Calculates also stiffness dSfdEf
% The sarcomere is embedded in Patch
% Theo Arts, Maastricht University, Oct 13, 2012

global P; % general time
Depolarization
Sarc= P.Patch;

%==== Input variables
t       = P.t;
Ef      = Sarc.Ef ;
tc      = repmat(t,[1,Sarc.n])-repmat(Sarc.Depolarization(1,:),[length(t),1]);
Lsi     = Sarc.Lsi;
C       = Sarc.C;

LenSeriesElement= Sarc.LenSeriesElement;
TR      = Sarc.TR     ; %TR atrial muscle > ventricular muscle
TD      = Sarc.TD     ;
TimeAct = Sarc.TimeAct; %time scale of contraction pulse
Ls0     = Sarc.Lsi0Act; %zero active stress sarcomere length
Ls0Pas  = Sarc.Ls0Pas ;
dLsPas  = Sarc.dLsPas ;
SfPas   = Sarc.SfPas  ;
CRest   = Sarc.CRest  ; %Resting C-value (Ca++  contractility)
LsRef   = Sarc.LsRef  ;
SfAct   = Sarc.SfAct  ;
vMax    = Sarc.vMax   ;

% series elasticity and sarcomere shortening velocity
Ls         = bsxfun(@times,exp(Ef),LsRef);
Sarc.Ls    = Ls;

%=== Active sarcomere
% constants related to timing are mainly based on experimental findings
L  = max(bsxfun(@rdivide,Lsi,Ls0)-1,0.0001) ; % normalized sarc length for active contraction
tA = bsxfun(@times,0.65+1.0570*L,TimeAct); % activation time lengthens with sarcomere length
tR = 0.55*TR.*TimeAct        ; % rise time
tD = 0.33*TD.*TimeAct        ; % decay time (default 0.22)
T  = bsxfun(@rdivide,tc,tR);
x  = min(8,max(0,T)); % normalized time during rise of activation
ft1= bsxfun(@rdivide, x.^3 .* exp(-x) .* (8-x).^2 * 0.020, tR);
%         rise of contraction, 'amount of Ca++ release'
%Integral T^n exp(-T) = Gamma(n+1) = n!
x= bsxfun(@rdivide,tc-tA,tD); % normalized time during decay of activation
tanhx= 0.5+0.5*sin( sign(x).*min(pi/2,abs(x)) ); %always>0
% Time confined approximation of 1/(1-e^x) function
FL= tanh(9.1204*L.^2); % regulates increase of contractility with Ls
Sarc.CDot= FL.*ft1 - bsxfun( @rdivide, C.*tanhx, tD); % 1st order rise and decay of [Ca++]
SfIso    = bsxfun(@times, C .* L,   1.51*SfAct) ;
SfRest   = bsxfun(@times, L, 1.51*CRest.*SfAct) ;

k1= 10; k2= 0.01; kk3= 2*(LsRef./dLsPas);
LfP   = bsxfun(@times,exp(Ef), LsRef./Ls0Pas);
yEcm  = LfP.^k1;
SfEcm = bsxfun(@times, yEcm-1, 0.0349*SfPas);% Ls=Ls0Pas, zero stress
y     = bsxfun(@power,LfP,kk3);
SfTit = bsxfun(@times, y-1, k2*SfAct);% titin is softer than ecm, proportional with SfAct
SfPasT= SfTit + SfEcm;
DSfPasDEf= bsxfun(@times,y, k2*SfAct.*kk3) + ...
    bsxfun(@times,yEcm, 0.0349*k1*SfPas);

%=== Stress/Ls and stiffness, collected for output/problem Buckle!!+++++
Sarc.SfEcm  = SfEcm;% passive stress ECM
Sarc.SfTit  = SfTit;% passive stress related to sarcomere: SfTit
LNormSe     = bsxfun(@rdivide,Ls-Lsi,LenSeriesElement);
Sarc.LsiDot = bsxfun(@times,LNormSe-1,vMax);
Sarc.Sf     = SfPasT + (SfIso+SfRest).*LNormSe - SfRest;
Sarc.DSfDEf = DSfPasDEf+ bsxfun(@rdivide,SfIso.*Ls,LenSeriesElement); % estimate of sarcomere stiffness

P.Patch=Sarc;
end

function Depolarization
% New time of depolarization
global P
Large    = 100; %used as infinity for timing

tDep     = P.Patch.Depolarization(1,:); %time of last depolarization
tNext    = P.Patch.Depolarization(2,:); %time of upcoming depolarization
DepPath  = P.Patch.DepPath; %paths of depolarization with delays

t        = P.t(end);
Rg       = tNext<t;
tDep(Rg) = tNext(Rg);
tNext(Rg)= t+Large;
tIn      = tDep(DepPath(1,:));
tOut1    = tIn + DepPath(3,:);
tOut     = tNext(DepPath(2,:));
tRefr    = tDep+P.Patch.TauRefrac;
tRefrac  = tRefr(DepPath(2,:));
Rg2      = (tOut1>tRefrac) & (tOut1<tOut);
tNext(DepPath(2,Rg2))= tOut1(Rg2);
if Rg2(P.Patch.iPace) % Cycle ends with activation of Pace-Patch
    P.tDot(end)=0;    % sign of  end of Cycle -> OdeCA
end
P.Patch.Depolarization=[tDep;tNext];
end

