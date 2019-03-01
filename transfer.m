clear;
filepattern = fullfile('E:\Data\HU', '*.bmp');
files = dir(filepattern);

for k=1:length(files)
    file = fullfile(files(k).folder, files(k).name);
    originPic=imread(file);
    numberOfpixel=350;
    GrayPic1=rgb2gray(int32(originPic));
    GrayPic2=imbinarize(GrayPic1,0.2);
    phi=mean(GrayPic1(:));

    GrayPic=GrayPic1-phi;

    FGrayPic=fft2(GrayPic);
    FGrayPic=fftshift(FGrayPic);
    MFGrayPic=abs(FGrayPic);
    MFGrayPic=MFGrayPic.^2;

    FlattenedData = MFGrayPic(:)'; 
    MappedFlattened = mapminmax(FlattenedData, 0, 255); 
    MappedData = reshape(MappedFlattened, size(MFGrayPic)); 

    sk=zeros(numberOfpixel*numberOfpixel,2);
    n=1;c=1;r=1;
    for i=1:1:numberOfpixel
        for j=1:1:numberOfpixel
            sk(n,1)=sqrt((i-(numberOfpixel/2))^2+(j-(numberOfpixel/2))^2);
            sk(n,2)=MappedData(i,j);
            n=n+1;
        end
    end
    func=sortrows(sk,1);
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

    figure;
    scatter(x,y,'k');
    axis([0 100 0 90]);
    set(gca,'FontName','Times New Roman','FontSize',18);
    box on;
    xlabel('$k$','Interpreter','latex');
    ylabel('S(k)','Rotation',0,'fontsize',18,'Position',[-20,45])
    %ylabel('$\hat \psi(k)$','Interpreter','latex','Rotation',0,'fontsize',18,'Position',[-20,45]);
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    saveas(gca,strcat('E:\Data\angular_averaged_', files(k).name));


    %figure;
    %imshow(originPic);

    %figure;
    %imshow(MappedData);
    I = MappedData;
    J = imadjust(I);
    figure, imshow(J)
    saveas(gca,strcat('E:\Data\2D_', files(k).name));
end