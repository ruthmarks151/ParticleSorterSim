function [] = ParticleSorterSim ()
%      y
%      |
%      Z--x
%     /
%[x y z]
simTime=2*10^-8;%also seconds
steps=500;
deltaT=simTime/steps;%seconds

t=0;


%B field is uniform
c=2.99*10^8;
bAMagnitude=[0,0,1];%units are in teslas
bAArea=[1,-0.17,-0.5;1.54,0.17,.5;];

bBMagnitude=[0,0,-1];%units are in teslas
bBArea=[3,-0.5,-0.5;3.54,-0.2,.5;];

pType=2;
%1 Electron
%2 Proton
%3 Pion
switch pType
    case 1
        chargeMassRatio=1.75*10^11;%in coulombs per kilogram
        vParticle=[.96*c,0,0];
        drawAs='b*';
    case 2
        chargeMassRatio=9.579*10^7;%in coulombs per kilogram
        vParticle=[.96*c,0,0];
        drawAs='r*';       
    case 3
        chargeMassRatio=6.439*10^8;%in coulombs per kilogram
        vParticle=[.96*c,0,0];
        drawAs='g*';
    otherwise
    
end

sParticle=[0,0,0];
hold on;
axis equal;
drawBField(bAArea,bAMagnitude);
drawBField(bBArea,bBMagnitude);

while(t<simTime)
    
    if inside(sParticle,bAArea)
        a=chargeMassRatio*cross(vParticle,bAMagnitude);
    elseif inside(sParticle,bBArea)
        a=chargeMassRatio*cross(vParticle,bBMagnitude);
    else
        a = [0,0,0];
    end
    sParticle=sParticle+vParticle*deltaT+(1/2)*a*deltaT*deltaT;
    vParticle=vParticle+a*deltaT;
    t=t+deltaT;
    plot3(sParticle(1),sParticle(2),sParticle(3),drawAs);
    drawnow;
end
end 

function [isInside] = inside(point,boundingArea)
    isInside=true;
    for i=1:3
        if point(i)>boundingArea(1,i)&&point(i)>boundingArea(2,i)
            isInside=false;
        end   
        if point(i)<boundingArea(1,i)&&point(i)<boundingArea(2,i)
            isInside=false;
        end  
    end
end

function [] = drawBField(bArea,bMagnitude)
    arrowsPerAxis=2;
    bMagnitude=bMagnitude*0.5/norm(bMagnitude);
    for x=min(bArea(:,1)):(max(bArea(:,1))-min(bArea(:,1)))/arrowsPerAxis:max(bArea(:,1))
        for y=min(bArea(:,2)):(max(bArea(:,2))-min(bArea(:,2)))/arrowsPerAxis:max(bArea(:,2))
            for z=min(bArea(:,3)):(max(bArea(:,3))-min(bArea(:,3)))/arrowsPerAxis:max(bArea(:,3))
                quiver3(x,y,z,bMagnitude(1),bMagnitude(2),bMagnitude(3),'black');
            end
        end
    end
end
