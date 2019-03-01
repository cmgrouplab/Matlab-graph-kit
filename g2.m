N=495; %Total number of particles
L=1; %The length of box
A=L*L; %The area of box
Rho=N/A; %Density of number of particles 

%Reading data from files
table1 = readtable('~/Documents/Results/Simulation2/1/1.txt'); %The coordinates of particle 1
table2 = readtable('~/Documents/Results/Simulation2/1/2.txt');
table3 = readtable('~/Documents/Results/Simulation2/1/3.txt');
table4 = readtable('~/Documents/Results/Simulation2/1/4.txt');
table5 = readtable('~/Documents/Results/Simulation2/1/5.txt');

%Create a new table for saving coordinates of all particles
table=[table1;table2;table3;table4;table5]; 

N_D=(N*(N-1))/2; %The number of distances (How many different distances)
D=zeros(N_D,1); %Create a matrix to record distances between any two particles

k=1;
for i=1:1:N   %Calculate the distance 
    for j=i+1:1:N
        tempx=abs(table.Var1(i)-table.Var1(j)); 
        tempy=abs(table.Var2(i)-table.Var2(j));
        if tempx>=L/2            %For boundary conditions
            tempx=L-tempx;
        end
        if tempy>=L/2
            tempy=L-tempy;
        end
        D(k)=sqrt(tempx*tempx+tempy*tempy); %Calculate distance between two particles
        k=k+1;
    end
end

m=ceil(((L/2)-0.01)/0.01)+1;
N_R=zeros(m,1); %The number of particles at R
x=zeros(1,m);

g2_R=zeros(1,m); %Create a matrix to save the value of g2

s=1;
for R=0.01:0.01:(L/2)  
    delta_R=0.001*R;
    for i=1:1:N_D
        if D(i)<=(R+(delta_R/2)) && D(i)>=(R-(delta_R/2))
            N_R(s)=N_R(s)+2;
        end
    end
    g2_R(s)=N_R(s)/(N*Rho*2*pi*R*delta_R);   %Calculate g2
    x(s)=R;
    s=s+1;
end
   
plot(x,g2_R,'bo:')
title('Simulation 2 Rho:0.55')
xlabel('R/L');
ylabel('g2(R)');
%saveas(gcf,'~/Desktop/Simulation2/s2-g2-16.jpg');

    
    
    
    