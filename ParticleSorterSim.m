function [] = ParticleSorterSim ()
%      y
%      |
%      Z--x
%     /
%[x y z]
%set up GPU array

DRAW_EACH_PARTICLE_INDIVIDUALLY = 1;

%Define details of time in the simulation
simTime=2*10^-8;%also seconds
steps=1000;
drawEvery = 10;
deltaT=simTime/steps;%seconds
t=0;
%define some handy values
c=3*10^8;
%Details about particles mass/charge ratio in coulombs per kilogram
electonCharge=-1.602*10^-19;
electronMass=9.109*10^-31;

protonCharge=1.602*10^-19;
protonMass=1.672*10^-27;

pionPositiveCharge = 1.602e-19;
pionNegativeCharge = -pionPositiveCharge;
pionChargedMass = 2.488064e-28;
pionNeutralCharge = 0;
pionNeutralMass = 2.406176e-28;

%Define All the B-Fields the areas are defined by two oppisite points
%The areas are in meters the field intensities are in Teslas
bAMagnitude=[0,0,1];
bAArea=[1,-0.17,-0.5;1.54,0.17,.5;];

bBMagnitude=[0,0,-1];
bBArea=[3,-0.5,-0.5;3.54,-0.2,.5;];

%particle source prefs
nParticles=100;
possibleParticles = {'proton','electron','pionPositive','pionNegative','pionNeutral'};
possibleParticlesWeights = [1,0,1,1,1];
velocitySpray=0.02*c;
%plot the b fields
hold on;
axis equal;
drawBField(bAArea,bAMagnitude);
drawBField(bBArea,bBMagnitude);
    

%Set up empty vectors
particleCount=0;
position=zeros(nParticles*1.1,3);
velocity=zeros(nParticles*1.1,3);
charge=zeros(nParticles*1.1);
mass=zeros(nParticles*1.1);
particleTypes = cell(1,nParticles * 1.1);

iterationNo = 0;
while(t<simTime)
    if random('unif',0,steps)< nParticles
        %initialize particles
        particleCount=particleCount+1;
        position(particleCount,:)=[0,0,0];
        velocity(particleCount,:)=[.96*c,random('unif',0,velocitySpray),random('unif',0,velocitySpray)];
        pType = datasample(possibleParticles,1, 'Weights', possibleParticlesWeights);
        particleTypes{particleCount} = pType{1};
        switch particleTypes{particleCount}
            case 'proton'
                charge(particleCount) = protonCharge;
                mass(particleCount) = protonMass;
            case 'electron'
                charge(particleCount) = electonCharge;
                mass(particleCount) = electronMass;
            case 'pionPositive'
                charge(particleCount) = pionPositiveCharge;
                mass(particleCount) = pionChargedMass;
            case 'pionNegative'
                charge(particleCount) = pionNegativeCharge;
                mass(particleCount) = pionChargedMass;
            case 'pionNeutral'
                charge(particleCount) = pionNeutralCharge;
                mass(particleCount) = pionNeutralMass;
        end
       
    end
    
    parfor id=1:particleCount
        if inside(position(id,:),bAArea)
            a=charge(id)/mass(id)*cross(velocity(id,:),bAMagnitude);
        elseif inside(position(id,:),bBArea)
            a=charge(id)/mass(id)*cross(velocity(id,:),bBMagnitude);
        else
            a = [0,0,0];
        end
        position(id,:)=position(id,:)+velocity(id,:)*deltaT+(1/2)*a*deltaT*deltaT;
        velocity(id,:)=velocity(id,:)+a*deltaT;
    end
    
    t=t+deltaT;
    if (mod(iterationNo, drawEvery) == 0)
        if DRAW_EACH_PARTICLE_INDIVIDUALLY
            for i=1:particleCount
                plot3(position(i,1), position(i,2), position(i,3),getDotType(particleTypes{i}));
            end
        else
            plot3(position(:,1),position(:,2),position(:,3),'bo');
        end
    end
    %particleCount
    drawnow;
    iterationNo = iterationNo + 1;
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
    bMagnitude=bMagnitude*0.25/norm(bMagnitude);
    for x=min(bArea(:,1)):(max(bArea(:,1))-min(bArea(:,1)))/arrowsPerAxis:max(bArea(:,1))
        for y=min(bArea(:,2)):(max(bArea(:,2))-min(bArea(:,2)))/arrowsPerAxis:max(bArea(:,2))
            for z=min(bArea(:,3)):(max(bArea(:,3))-min(bArea(:,3)))/arrowsPerAxis:max(bArea(:,3))
                quiver3(x,y,z,bMagnitude(1),bMagnitude(2),bMagnitude(3),'black');
            end
        end
    end
end

function [dotType] = getDotType(particleType)
    dotType = 'blacko';
    switch particleType
        case 'proton'
            dotType = 'ro';
        case 'electron'
            dotType = 'bo';
        case 'pionPositive'
            dotType = 'r*';
        case 'pionNegative'
            dotType = 'b*';
        case 'pionNeutral'
            dotType = 'black*';
    end
end