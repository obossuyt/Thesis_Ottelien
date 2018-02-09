function [deltatot,V]= funPark_atSpecificLocationHeight(X,Y,windfarm,windturbine,wake,k,options,H)
N=length(windfarm.WTlocx);
D=windturbine.D;
Hhub=windturbine.Hhub;

% Define in which wakes every point is located
for i=1:N
    [wake(i).in,wake(i).on]=inpolygon(X,Y,wake(i).xWake,wake(i).yWake);
end

delta=0; %initialize wake deficit to 0
%         scatter(x(idx),y(idy));

dh=Hhub-H;
for i=1:N
    if wake(i).in==0
        delta=delta; %no extra wake deficit
    else
        % Compute distance perpendicular to rotor plane
        curvexy=[wake(i).xCenterline',wake(i).yCenterline'];
        mapxy=[X,Y];
        [xy,~,~] = distance2curve(curvexy,mapxy,'linear');
        %                                                  scatter(xy(1),xy(2),'*');
        XCwake=xy(1);
        YCwake=xy(2);
        
        %Compute contribution of each wake
        if i==1 || options.ParkModel==2
            term1=1-sqrt(1-wake(i).CT);
        else
            term1= 1- wake(i).V/wake(i).U*sqrt(1-wake(i).CT);
        end
        dist=sqrt((XCwake-windfarm.WTlocx(i))^2+(YCwake-windfarm.WTlocy(i))^2);
        Rwake=1/2*(D+2*k*dist);
        c=2*(Rwake^2-dh^2)^(1/2);
        
     if real(c)>0  && sqrt( (X-XCwake)^2+ (Y-YCwake)^2) < c/2

        
        
        
        term2=(D/(2*Rwake))^2;
        
        if options.AreaAveraging==1
            Aoverlap=AreaOverlap([Y,windfarm.WTlocy(i)], [H,Hhub], [D/2,c/2]);
            term3=Aoverlap/wake(i).A;
        else
            term3=1;
        end
        res=term1*term2*term3;
        
       if k*dist>2*Hhub-D && options.WakeReflection % Add deficit due to wake reflection as well
            if strcmp(options.SPmethod,'lin')
                delta=delta+2*res;
            elseif strcmp(options.SPmethod,'quadr')
                delta=delta+2*res^2;
            elseif strcmp(options.SPmethod,'max')
                delta=max(delta,res);
            end
        else
            if strcmp(options.SPmethod,'lin')
                delta=delta+res;
            elseif strcmp(options.SPmethod,'quadr')
                delta=delta+res^2;
            elseif strcmp(options.SPmethod,'max')
                delta=max(delta,res);
                
            end
        end
     else
         delta=delta;
     end
     
        
    end
end
if strcmp(options.SPmethod,'quadr')
    deltatot=sqrt(delta);
else
    deltatot=delta;
end

V=wake(i).U*(1-deltatot);