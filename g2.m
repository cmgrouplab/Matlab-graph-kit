clear;
N=274; %Total number of particles
L=1; %The length of box
A=L*L; %It is the area of box
Rho=N/A; %Density of number of particles 


%table1 = readtable('~/Desktop/5.1/0.5/1.txt'); %Reading data from files
%table2 = readtable('~/Desktop/5.1/0.5/2.txt');

%table=[table1;table2];
table = readtable('centroidPositions_2D_amorphous_SiO2_particle.txt');
table.Var1 = table.Var1 * L;
table.Var2 = table.Var2 * L;

N_D=(N*(N-1))/2; %The number of distances
D=zeros(N_D,1); %Create a matrix to record distances between any two particles

k=1;
for i=1:1:N   %Calculate the distance 
    for j=i+1:1:N
        tempx=abs(table.Var1(i)-table.Var1(j));
        tempy=abs(table.Var2(i)-table.Var2(j));
        if tempx>=L/2
            tempx=L-tempx;
        end
        if tempy>=L/2
            tempy=L-tempy;
        end
        D(k)=sqrt(tempx*tempx+tempy*tempy); %Calculate distance between two particles
        k=k+1;
    end
end

m=ceil(((L/2)-0.01)/0.015)+1;
N_R=zeros(m,1); %The number of particles at R
x=zeros(1,m);

g2_R=zeros(1,m);

s=1;
for R=0.01:0.015:(L/2)  
    delta_R=0.05*R;
    for i=1:1:N_D
        if D(i)<=(R+(delta_R/2)) && D(i)>=(R-(delta_R/2))
            N_R(s)=N_R(s)+2;
        end
    end
    g2_R(s)=N_R(s)/(N*Rho*2*pi*R*delta_R);
    x(s)=R;
    s=s+1;
end

x1 = x/L;
g2Data = [x1;g2_R];
g2Data = g2Data';
plot(x1(1:33),g2_R(1:33),'ko:')
title('g2 function for 2D amorphous SiO2 particle')
xlabel('R/L');
ylabel('g2(R)');
%axis([0 0.5 0 3]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
 box on;
 saveas(gca,'g2_function_for_2D_amorphous_SiO2_particle.png');
 save g2_function_for_2D_amorphous_SiO2_particle.txt -ascii g2Data



    
    
    
    