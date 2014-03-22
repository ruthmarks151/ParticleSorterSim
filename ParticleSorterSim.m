%      y
%      |
%      O--x
%     /
%    Z
%[x y z]
simTime=1*10^-8;%also seconds
steps=1000
deltaT=simTime/steps;%seconds

t=0;


%B field is uniform
c=2.99*10^8;
b=[0,0,0*10^-3];%units are in teslas
%Electric field
e=[0,1,0];
magLength=0.56;

pType=1;
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

sParticle=[0,0.5,0];
hold on;

while(t<simTime)
    if sParticle(1)<magLength
        a=chargeMassRatio*cross(vParticle,b)+chargeMassRatio*e;
    else
        a=[0,0,0];
    end
    sParticle=sParticle+vParticle*deltaT+(1/2)*a*deltaT*deltaT;
    vParticle=vParticle+a*deltaT;
    t=t+deltaT;
    plot3(sParticle(1),sParticle(2),sParticle(3),drawAs);
    drawnow;
end
