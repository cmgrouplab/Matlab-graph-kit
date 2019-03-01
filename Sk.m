Length=1; %length of box
Number=1000; %number of particles;
c=1;r=1;

n=1:20;
kx(1,:)=n*((2*pi)/Length);
[kx,ky,kz] = ndgrid(kx,kx,kx);
k=[kx(:),ky(:),kz(:)];

A=importdata('/Users/yuzheng/Documents/Results/LSD.HS.mono.3D.N=1000.JamT1.001.lsdpg.dat',' ',6);
positions=A.data;
positions=positions';

product=k*positions;
exp1=exp(1i*product);

sum1=sum(exp1,2);
kmod=abs(sum1);
sk=(kmod.^2)/Number;

k1=k.^2;    
k2=sum(k1,2);
ksqrt=k2.^(1/2);

func=[ksqrt,sk];
func=sortrows(func,1);
B=tabulate(func(:,1));
[row,col]=size(func);
[rowb,colb]=size(B);
final=zeros(rowb,2);
while c<row+1
    d=B(r,2);
    final(r,:)=(sum(func(c:(c+d-1),:),1))/d;
    c=c+d;
    r=r+1;
end
finals=final';

x(1,:)=finals(1,:);
y(1,:)=finals(2,:);

scatter(x,y,'k');
title('Sk.JamT1.001');
xlabel('k');
ylabel('S(k)');
%saveas(gcf,'~/Desktop/Sk.JamT1.png');
    