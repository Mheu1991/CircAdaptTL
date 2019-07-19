function CircDisplay(Legends)
% function CircDisplay(Legends)
% Displays hemodynamics, based on structure P
% If Legends==TRUE, legends are added
% Theo Arts, Maastricht University, Feb 29, 2016

global P

SVarDot(0,P.SVar',[]);
if ~exist('Legends')
    Legends=0;
end

% scaling of graphics
q0          = P.General.q0;
p0          = P.General.p0;
tCycle      = P.General.tCycle;

%=== scaling reference values
qSc= Rnd(q0);
VSc= Rnd(q0*tCycle/10);
pSc= Rnd(0.1*p0);

t= P.t-P.t(1);
OFFSET= -20;

% find arterial heart valves with distal nodes and valves proximal to atria
% resulting strings are used to identify nodes and valves
nValve=P.Valve.n;
for i=1:nValve
    if regexpi(P.Valve.Name{i},'Lv')==1
        LvAo= P.Valve.Name{i};
        Ao  = P.Node.Name(P.Valve.iNodeDist(i));
    end;
    if regexpi(P.Valve.Name{i},'Rv')==1
        RvPu= P.Valve.Name{i};
        Pu  =P.Node.Name(P.Valve.iNodeDist(i));
    end;
    if regexpi(P.Valve.Name{i},'La')>1
        VeLa= P.Valve.Name{i};
    end;
    if regexpi(P.Valve.Name{i},'Ra')>1
        VeRa= P.Valve.Name{i};
    end;
end

%==== Hemodynamics
figure(1);
p1=Get('Node','p',{Ao,'Lv',Pu,'Rv'})/pSc;
V1=Get('Cavity','V',{'La','Lv','Ra','Rv'})/VSc;
q1=Get('Valve','q',{LvAo,'LaLv',VeLa,RvPu,'RaRv',VeRa})/qSc;
p1(:,[3,4])= p1(:,[3,4])+OFFSET;% pressures
V1(:,[3,4])= V1(:,[3,4])+OFFSET;% volumes
q1(:,4:end)= q1(:,4:end)+OFFSET;% flows
A=[p1,V1,q1];
Maxx=max(t);
Minx=min(t);
Maxy=max(A(:));
Miny=min(A(:));
subplot(2,4,[3,4,7,8]); plot(t,A);
axis([Minx,Maxx,Miny,Maxy]);
title(['Units: ',num2str(qSc*1e6),' ml/s; ',num2str(pSc/1e3),' kPa; ',...
    num2str(VSc*1e6),' ml']);

% p-V loops
Calp= 0.001; CalV= 1e6;
VT= CalV*Get('Cavity','V',{'La','Ra','Lv','Rv','Rv'});
pT= Calp*Get('Cavity','p',{'La','Ra','Lv','Rv','Rv'});
Maxx=max(VT(:));
Minx=min(VT(:));
Maxy=max(pT(:));
Miny=min(pT(:));
subplot(2,4,1); plot(VT,pT);
axis([0,Maxx,Miny,Maxy]);
title('p(V)')
if Legends
    legend({'La','Ra','Lv','Rv'})
end;


%Sf-Ef loops
EfT = Get('Patch','Ef','All');
SfT = Get('Patch','Sf','All')/1e3;
Maxx=max(EfT(:));
Minx=min(EfT(:));
Maxy=max(SfT(:));
Miny=min(SfT(:));
subplot(2,4,2); plot(EfT,SfT);
axis([Minx,Maxx,Miny,Maxy]);
title('stress[kPa](strain)')
if Legends
    legend(P.Patch.Name)
end;

% clipped low venous, atrial and ventricular pressures
A= [Get('Cavity','p',{'La','Ra','Lv','Rv','Rv'}),...
    Get('Bag'   ,'p',{'Peri'             })]/1e3;
Maxx=max(t);
Minx=min(t);
Maxy=max(max(A(:,[1,2])));
Miny=min(A(:));
subplot(2,4,5:6); plot(t,A);
axis([Minx,Maxx,Miny,Maxy]);
title('p(t) [kPa]')
if Legends
    legend({'La','Ra','Lv','Rv','Peri'})
end;

figure(1);
end


%========== AUXILARY FUNCTIONS =============

function X=Rnd(x)
X1= 10^round(log(x)/log(10));
X2=  2^round(log(x/X1)/log(2));
X=X1*X2;
end

