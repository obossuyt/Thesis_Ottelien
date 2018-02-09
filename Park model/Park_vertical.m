clear all
close all

%% DEFINE METHODS AND LOAD CASE
options.LoadCase=1; % 1, 2 or 3
options.SPmethod='max'; % Superposition method: 'lin'= linear SP, 'max'=maximum deficit, 'quadr'=quadratic SP
options.AreaAveraging=false;

options.ParkModel=1; %1: Park 1 model (Correctionfactor in front of sqrt(1-CT)); 2: Park 2 model
options.WakeReflection=false; % True: include wake reflection; False: no wake reflection

options.MakePlots=true;

k=0.04;             % Wake decay coefficient

%% DATA
abl0vec = [236.3 243.1 236.7]; % From Nicolai, ASSUMED TRUE INFLOW DIRECTION AT HUB HEIGHT
U0vec = [11.8 9.0 9.5]; % Free stream velocity

theta0=abl0vec(options.LoadCase);
U0=U0vec(options.LoadCase);

D=154;
Hhub=106;

%% SET locations

% Downstream distance
X1=20*D;

% Cross-stream distance
Y2=0*D;

%% WINDTURBINE: wind turbine data

WTdata=xlsread('TurbineData.xlsx','Turbine2rescaled');
data=load('TurbineData.mat');
WTdata=data.WTdata;
WSvec=WTdata(:,1);
[windturbine.WSvec,idx]=unique(WSvec);
windturbine.CPvec=WTdata(idx,2);
windturbine.CTvec=WTdata(idx,3);
windturbine.D=D;               % Rotor diameter
windturbine.Hhub=Hhub;

%% WINDFARM: Wind turbine coordinates and data
% Westermost Rough Wind Farm
dataCo=load('WTcoord.mat');
dataCo=dataCo.WTcoord;
WTx=[dataCo([1 8 15],1)]; % [1 8 15]
WTy=[dataCo([1 8 15],2)];

% WTx=[dataCo([1 8 15 20],1)]; % [1 8 15]
% WTy=[dataCo([1 8 15 20],2)];

WTx=WTx-WTx(1); % move to x=0 for WT A01
WTy=WTy-WTy(1); % move to y=0 for WT A01

% Rotate wind turbine so wind is coming from the west
thetadeg= -(270-theta0);
windfarm.WTlocx = WTx*cosd(thetadeg) - WTy*sind(thetadeg);
windfarm.WTlocy = WTx*sind(thetadeg)+ WTy*cosd(thetadeg);

% Show as wind farm is located in reality: Adapt theta!
% windfarm.WTlocx=WTx;
% windfarm.WTlocy=WTy;
% windfarm.theta=theta0*ones(1,length(windfarm.WTlocx)); %Wind direction in windDir coord (N=0°, E=90°, S=180°, W=270°)

% Sort from upstream to downstream
A=[windfarm.WTlocx,windfarm.WTlocy];
B=sortrows(A,1);
windfarm.WTlocx=B(:,1);
windfarm.WTlocy=B(:,2);

N=length(windfarm.WTlocx);
% if (N ~= length(windfarm.WTlocy))
%     error('Mistake in defining wind farm! WTlocx and WTlocy should have the same length.')
% end

% Wind direction in windDir coord (N=0°, E=90°, S=180°, W=270°)
windfarm.theta=270*ones(1,length(windfarm.WTlocx));
% yawerror=load('YawError.mat');
% yawerror=yawerror.ell(options.LoadCase).yaw;
% windfarm.theta=windfarm.theta-yawerror';


%% Grid setup
% Grid around every wind turbine
dx=0.1*D;
dy=0.05*D;
dz=0.02*D;

x=-2*D+min(windfarm.WTlocx):dx:(10*D+max(windfarm.WTlocx));     % x spacing grid
Dmax=round((D+2*k*max(x))/D)+1;
y=-2*D+min(windfarm.WTlocy):dy:max(windfarm.WTlocy)+2*D;   % y spacing grid
Hmax=round(Hhub+1/2*(D+2*k*max(x)))/D;
z=0:dz:1.1*Hmax*D;

%% Calculation of individual wakes
wake = funPark_atWindTurbines(windfarm,windturbine,U0,k,x,options);

%% Cross-stream section at particular downstream position

[Y1,Z1]=meshgrid(y,z);

% Define downstream distance
% X1=14.77*D;

T=tic;
% Loop over all wanted locations
for idy=1:length(y)
    for idh=1:length(z)
        [deltatot1(idh,idy),V1(idh,idy)]=funPark_atSpecificLocationHeight(X1,y(idy),windfarm,windturbine,wake,k,options,z(idh));
    end
end
t=toc(T);

%% Downwind cross-section at particular cross-stream position

% Define cross-stream position
% Y2=60;

[X2,Z2]=meshgrid(x,z);
T=tic;
% Loop over all wanted locations
for idx=1:length(x)
    for idh=1:length(z)
        [deltatot2(idh,idx),V2(idh,idx)]=funPark_atSpecificLocationHeight(x(idx),Y2,windfarm,windturbine,wake,k,options,z(idh));
    end
end
t=toc(T);

%% Plots

if options.MakePlots
    
    %% At specific downstream distance
     figure
    hold on
    [C,ch]=contourf(Y1/D,Z1/D,deltatot1,50);
    set(ch,'LineColor','none')
    scatter(windfarm.WTlocy/D,Hhub*ones(N,1)/D,'k.')
    c=colorbar;
    c.Label.String='\delta U';
    xlabel('y/D')
    ylabel('z/D')
    axis equal
    grid on
     title(sprintf('x=%1.2fD',X1/D))
    
    figure
    hold on
    [C,ch]=contourf(Y1/D,Z1/D,V1,50);
    set(ch,'LineColor','none')
    scatter(windfarm.WTlocy/D,Hhub*ones(N,1)/D,'k.')
    c=colorbar;
    c.Label.String='U [m/s]';
    colormap jet
    xlabel('y/D')
    ylabel('z/D')
    axis equal
    grid on
    title(sprintf('x=%1.2fD',X1/D))
%     centers = [windfarm.WTlocy Hhub*ones(3,1)]/D;
%  radii = D/2*ones(3,1)/D;
% % % Display the circles.
%  viscircles(centers,radii,'color','k');
    
%%     At specific cross-stream position
   
    figure
    hold on
    [C,ch]=contourf(X2/D,Z2/D,deltatot2,50);
    set(ch,'LineColor','none')
    scatter(windfarm.WTlocy/D,Hhub*ones(N,1)/D,'k')
    c=colorbar;
    c.Label.String='\delta U';
    xlabel('x/D')
    ylabel('z/D')
    axis equal
    grid on
     title(sprintf('y=%1.2fD',Y2/D))
    
    figure
    hold on
    [C,ch]=contourf(X2/D,Z2/D,V2,50);
    set(ch,'LineColor','none')
    scatter(windfarm.WTlocx/D,Hhub*ones(N,1)/D,'k')
    c=colorbar;
    c.Label.String='U [m/s]';
    colormap jet
    xlabel('x/D')
    ylabel('z/D')
    axis equal
    grid on
     title(sprintf('y=%1.2fD',Y2/D))

end


%% Save figures
% filename=strcat(sprintf('Figures/WMR_theta%1.0f_',theta0),SPmethod,'_rotated');
% saveas(gcf,strcat(filename,'.png'))
% PlotToFileColorPDF(gcf,filename,10,10)

