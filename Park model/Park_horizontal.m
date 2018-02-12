clear all
close all

%% DEFINE METHODS AND LOAD CASE
options.LoadCase=1; % 1, 2 or 3
options.SPmethod='quadr'; % Superposition method: 'lin'= linear SP, 'max'=maximum deficit, 'quadr'=quadratic SP
options.AreaAveraging=false;

options.ParkModel=1; %1: Park 1 model (Correctionfactor in front of sqrt(1-CT)); 2: Park 2 model
options.WakeReflection=true; % True: include wake reflection; False: no wake reflection
options.WTdata=2; %1: LES simulations; 2: WMR turbine

options.MakePlots=true;

k=0.04;             % Wake decay coefficient

%% SET height for plot

H=106;

%% DATA
abl0vec = [236.3 243.1 236.7]; % From Nicolai, ASSUMED TRUE INFLOW DIRECTION AT HUB HEIGHT
U0vec = [11.8 9.0 9.5]; % Free stream velocity

theta0=abl0vec(options.LoadCase);
U0=U0vec(options.LoadCase);

%% WINDTURBINE: wind turbine data

D=154;
Hhub=106;

%WTdata=xlsread('TurbineData.xlsx','Turbine2rescaled');
if  options.WTdata==1
data=load('TurbineData.mat');
WTdata=data.WTdata;
elseif options.WTdata==2

data=load('WMRturbine.mat');
WTdata=data.WMRdata;
end


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
WTx=[dataCo([1 8 15 20],1)]; % [1 8 15]
WTy=[dataCo([1 8 15 20],2)];

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
dy=0.02*D;

x=-2*D+min(windfarm.WTlocx):dx:(10*D+max(windfarm.WTlocx));     % x spacing grid
Dmax=round((D+2*k*max(x))/D)+1;
y=-Dmax*D+min(windfarm.WTlocy):dy:max(windfarm.WTlocy)+Dmax*D;   % y spacing grid

[X,Y]=meshgrid(x,y);

%% Calculation of individual wakes
wake = funPark_atWindTurbines(windfarm,windturbine,U0,k,x,options);

%% Define wind speed at specific locations

T=tic;
% Loop over all wanted locations
for idx=1:length(x)
    for idy=1:length(y)
        [deltatot(idy,idx),V(idy,idx)]=funPark_atSpecificLocationHeight(x(idx),y(idy),windfarm,windturbine,wake,k,options,H);
    end
end
t=toc(T);

%% Plots

if options.MakePlots
    figure
    hold on
    [C,h]=contourf(X/D,Y/D,deltatot,50);
    set(h,'LineColor','none')
    scatter(windfarm.WTlocx/D,windfarm.WTlocy/D,'k','filled')
    c=colorbar;
    c.Label.String='\delta U';
    xlabel('x/D')
    ylabel('y/D')
    axis equal
    grid on
    title(sprintf('z=%1.2f m',H))
    
    figure
    hold on
    [C,h]=contourf(X/D,Y/D,V,50);
    set(h,'LineColor','none')
    scatter(windfarm.WTlocx/D,windfarm.WTlocy/D,'k','filled')
    plot(windfarm.WTlocx(1)/D*ones(1,2),windfarm.WTlocy/D+[D/2 -D/2]/D,'k','linewidth',2)
    plot(windfarm.WTlocx(2)/D*ones(1,2),windfarm.WTlocy/D+[D/2 -D/2]/D,'k','linewidth',2)
    plot(windfarm.WTlocx(3)/D*ones(1,2),windfarm.WTlocy/D+[D/2 -D/2]/D,'k','linewidth',2)
    c=colorbar;
    c.Label.String='U [m/s]';
    colormap jet
    xlabel('x/D')
    ylabel('y/D')
    axis equal
    grid on
    title(sprintf('z=%1.2f m',H))

    %  xlim([-2+min(WTlocx/D),max(WTlocx/D)+10])
    %  ylim([-2+min(WTlocy/D),max(WTlocy/D)+5])
end

%% Save figures
%  filename=strcat(sprintf('Figures/WMR_theta%1.0f_z%1.1fm_',theta0,H),options.SPmethod)
%  saveas(gcf,strcat(filename,'.png'))
% PlotToFileColorPDF(gcf,filename,10,10)

