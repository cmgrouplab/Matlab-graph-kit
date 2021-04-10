clear;
Length=1; %length of box
Number=79; %number of particles;
c=1;r=1;

n=1:1:100;
kx(1,:)=n*((2*pi)/Length);
[kx,ky] = ndgrid(kx,kx);
k=[kx(:),ky(:)];
%k = [kx(:)];
A=importdata('1.txt');
positions=A;
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
x1 = x(1:200);
y1 = y(1:200);

skData = [x1;y1];
skData = skData';

scatter(x1,y1,600,'.','k');
title('Sk for 2D Material');
xlabel('k');
ylabel('S(k)');
box on;
axis([0 150 0 10])
saveas(gcf,'sk_2D_Material.png');
 save sk_for_2D_Material.txt -ascii skData
    