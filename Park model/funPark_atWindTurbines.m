function wake = funPark_atWindTurbines(windfarm,windturbine,U0,k,x,options)


% SPmethod=options.SPmethod;
% ParkModel=options.ParkModel;

WTlocx=windfarm.WTlocx;
WTlocy=windfarm.WTlocy;
theta=windfarm.theta;
D=windturbine.D;
Hhub=windturbine.Hhub;
WSvec=windturbine.WSvec;
CTvec=windturbine.CTvec;




% Initialization for first wind turbine
wake(1).U=U0;
wake(1).V=U0;
wake(1).CT=interp1(WSvec,CTvec,wake(1).U);
wake(1).delta=0;

alpha=mod((270-theta),360); %Conversion to conventional degree formulation

figure
hold on
N=length(WTlocx);
for i=1:N
    
    % Rotor diameter
    wake(i).A=pi*D^2/4;
    
    % Wake diameter at end of wake
    Rmax=1/2*(D+2*k*(x(end)-WTlocx(i)));
    
    %Wake polygon and rotated wake polygon
    xW=[WTlocx(i), WTlocx(i), x(end), x(end), WTlocx(i)]-WTlocx(i);
    yW=[WTlocy(i)-D/2,WTlocy(i)+ D/2, WTlocy(i)+Rmax,WTlocy(i)-Rmax,WTlocy(i)-D/2]-WTlocy(i);
    rotM = [cosd(alpha(i)),-sind(alpha(i));sind(alpha(i)),cosd(alpha(i))];
    rData = rotM*[xW;yW]; % compute 2-by-N array of rotated data
    
    % Rotated wake centerline
    rData2=rotM*[linspace(0,x(end)-WTlocx(i),length(x));...
        linspace(0,0,length(x))];
    wake(i).xCenterline=rData2(1,:)+WTlocx(i);
    wake(i).yCenterline=rData2(2,:)+WTlocy(i);
    
    % Polygon which defines wake area of each wind turbine
    wake(i).xWake = rData(1,:)+WTlocx(i)+0.0001; % extract x-coordinates of rotated data
    wake(i).yWake = rData(2,:)+WTlocy(i); % extract y-coordinates of rotated data
    
    plot(wake(i).xWake,wake(i).yWake);
    hold on
    plot(wake(i).xCenterline,wake(i).yCenterline,'--');
    
end
scatter(WTlocx,WTlocy,'k')
xlabel('x')
ylabel('y')
grid on

%% Define wind speeds at wind turbine -- TESTED WITH SIMPLE TEST CASE, SHOULD BE OK
for i=2:N
    wake(i).deltatot=0;
    for j=1:i-1 %loop over upstream wind turbines
        %Normal distance to upstream wind turbine
        xdist=WTlocx(i)-WTlocx(j);
        ydist=WTlocy(i)-WTlocy(j);
        d=sqrt(xdist^2+ydist^2);
        beta=atand(ydist/xdist);
        gamma=beta-alpha(j);
        wake(i).X(j)=d*cosd(gamma);
        
        
        temp=D+2*k*wake(i).X(j); %wake diameter at downstream wind turbine location
        wake(i).Aoverlap(j)= AreaOverlap( [WTlocy(j) WTlocy(j)-sqrt(d^2-wake(i).X(j)^2)], [Hhub Hhub], [temp/2 , D/2]);
        if options.ParkModel==1
            term1=1-wake(j).V/wake(j).U*sqrt(1-wake(j).CT);
        else %ParkModel==2
            term1=1-sqrt(1-wake(j).CT);
        end
        term2=(D/(D+2*k*wake(i).X(j)))^2;
        term3=wake(i).Aoverlap(j)/wake(i).A;
        wake(i).delta(j)=term1*term2*term3;
        
        if k*wake(i).X(j)>2*Hhub-D && options.WakeReflection % Add deficit due to wake reflection as well
            if strcmp(options.SPmethod,'lin')
                wake(i).deltatot=wake(i).deltatot+2*wake(i).delta(j);
            elseif strcmp(options.SPmethod,'quadr')
                wake(i).deltatot=wake(i).deltatot+2*wake(i).delta(j)^2;
            elseif strcmp(options.SPmethod,'max')
                wake(i).deltatot=max(wake(i).deltatot, wake(i).delta(j));     
            end
        else
            if strcmp(options.SPmethod,'lin')
                wake(i).deltatot=wake(i).deltatot+wake(i).delta(j);
            elseif strcmp(options.SPmethod,'quadr')
                wake(i).deltatot=wake(i).deltatot+wake(i).delta(j)^2;
            elseif strcmp(options.SPmethod,'max')
                wake(i).deltatot=max(wake(i).deltatot, wake(i).delta(j));     
            end
        end
        
        wake(i).Vtemp(j)=wake(i-1).U * (1-wake(i).delta(j));
    end
    
    if strcmp(options.SPmethod,'quadr')
        wake(i).deltatot=sqrt(wake(i).deltatot);
    end
    wake(i).U=wake(i).Vtemp(i-1)/(1-wake(i).delta(i-1));
    wake(i).V=wake(i).U*(1-wake(i).deltatot);
    wake(i).CT=interp1(WSvec,CTvec,wake(i).V);
end