function Aoverlap = AreaOverlap(Y,Z,R)

% only vertical and/or horizontal distance matters, distance in streamwise
% direction does not.

d=sqrt( (Y(1)-Y(2))^2+(Z(1)-Z(2))^2);
r1=R(1);
r2=R(2);

C1C2 = d;
% figure
% centers = [Y' Z'];
% radii = [r1 r2];
% % Display the circles.
% viscircles(centers,radii);
% hold off
% axis square

%%

if C1C2== r1+r2
    %disp('Wakes touch to each other but do not overlap')
    Aoverlap=0;
elseif C1C2 > r1+r2
    %disp('Wakes do not overlap')
    Aoverlap=0;
elseif C1C2 < r1+r2
    
%     figure
% centers = [Y' Z'];
% radii = [r1 r2];
% % Display the circles.
% viscircles(centers,radii);
% hold off
% axis square

    %disp('Wakes overlap')
    if d==0
        Aoverlap = min(pi*r1^2,pi*r2^2);
    else       
        term1=r1^2*acos((d^2+r1^2-r2^2)/(2*d*r1));
        term2=r2^2*acos((d^2+r2^2-r1^2)/(2*d*r2));
        term3=1/2*sqrt( (-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2));
        Aoverlap=term1+term2-term3;
        if imag(Aoverlap)~=0
            Aoverlap = min(pi*r1^2,pi*r2^2);
        end
    end
end
