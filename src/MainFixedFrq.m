set(0,'DefaultFigureWindowStyle','docked')

% Cleaning
clear all
close all
clc


addpath('FGen')
%% Primary tests

kr = 0:0.0001:10;
figure,hold on,
for p = 0:100
    plot(kr,abs(besselj(p,kr)));
end



%% Parameters

global DiscreteDiskCentered Feta HEST List_k % discretization of the centered artificial disk

% Choose Geometry between the available geometries:
Rin = 0.01; %0.1; % radius of artificial background
% Folder = 'FixedFrqData/1Crack/';%'Geometry_6Cracks_k14.mat';%'Geometry_2Cracks_k14.mat';
% load([Folder 'FFRin01'  '.mat'],'Ax','Ay','Bx','By','FF','List_k');
% num2str(round(10*Rin))

for iterNbCracks = 2:8
FileInput = ['MultiFrqData/Monotonicity/FF' num2str(iterNbCracks) '.mat'];
load(FileInput,'Ax','Ay','Bx','By','List_k',...
'BXmin','BXmax','BYmin','BYmax','FF');

%%

[nbcapteur,~,Nk] = size(FF);
% Choose resolution of the scan
% 
% BXmin = -.55; % range of the centers coordinates
% BXmax = .45;
% BYmin = -.5;
% BYmax = .5;
DB = Rin; % distance between each centers
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);

% Choose artificial background boundary condition:
% BC1 u + BC2 d_n u = 0.
BC = [1 0];

% plot and/or save computed indicators?
PlotDisplayFlag = 1; 
DataSaveFlag = 1;
PlotSaveFlag = 0;

% Choose noise level of farfield in the next section (eta).
% 
% % Control of the scan parameters
% figControle = figure;%('visible','off');
% hold on, 
% plot([Ax;Bx],[Ay;By],'red');
% set(gca,'Ydir','normal')
% axis square
% colormap('summer');
% title('Resolution of the image');
% 
% 
% % visualisation of the artificial backgrounds (when possible):
% if Nt<1500
% poly = DiskSwipe(BX,BY,Rin,30);
% figure,hold on,
% plot(poly,'FaceAlpha',0.1,'EdgeAlpha',0);
% plot([Ax;Bx],[Ay;By],'red','LineWidth',1);
% axis square
% title('The artificial backgrounds');
% end

%% Add of noise on the data
clc
eta=.02;
Feta = zeros(nbcapteur,nbcapteur,Nk);
Fart = zeros(nbcapteur,nbcapteur,Nk);
HEST = zeros(1,Nk);

for i=1:Nk

     
    % Bruit uniforme
    % noise=1+eta*(rand(nbcapteur)*2-1)+sqrt(-1)*eta*(rand(nbcapteur)*2-1);

    % Bruit Gaussien
     noise=1+eta*(randn(size(FF(:,:,i))))+sqrt(-1)*eta*(randn(size(FF(:,:,i))));

    Feta(:,:,i)=FF(:,:,i).*noise;

    HEST(i)=norm(Feta(:,:,i)-FF(:,:,i));
    if HEST(i) < 1.e-02
        HEST(i)= 1.e-02;                       
    end
    
    nF=norm(FF(:,:,i));
    disp(['Bruit en pourcentage = ',num2str(HEST(i)/nF*100)])
end


%% Indicators computation



% I, J are two different indicators - their definition are provided in the 
% paper. The number 1 or 2 indicates the method % used to compute I or J 
% (1=Tikhonov, 2=GLSM_FSharp).
I1 = zeros([size(BX) Nk]); I2 = zeros([size(BX) Nk]); J1 = zeros([size(BX) Nk]); J2 = zeros([size(BX) Nk]);

I1Temp = zeros(size(BX));
I2Temp = zeros(size(BX));
J1Temp = zeros(size(BX));
J2Temp = zeros(size(BX));

% Sampling points on the centered artificial background
Nbz = 20; % Nb Approx de points de discretisation interieurs du disque artificiel
[X,Y]=meshgrid(linspace(-Rin,Rin,sqrt(4./pi*Nbz)));
Bool = X.^2+Y.^2 < (Rin).^2;
DiscreteDiskCentered = transpose([X(Bool),Y(Bool)]);
NbPoints = length(DiscreteDiskCentered) % Nb de points de discretisation interieurs du disque artificiel
% figure,plot(X(Bool),Y(Bool),'red+'), 
% hold on, plot(Rin*cos(0:0.1:7),Rin*sin(0:0.1:7),'black'),
% axis square;
% title(['Discretisation interieur du disque artificiel centre de rayon ' num2str(Rin) ', nbPtsInterieur: ' num2str(NbPoints)]);
                                                                                                             



for iter_k = Nk %1:Nk % Boucles sur k 
    tic
        disp(['k= ' num2str(iter_k) '/' num2str(Nk)]);
        k = List_k(iter_k);
        
        
%Etalonnage 
tCalibration = [-10,-10];
% FDiskCentered = FGen(k,Rin,[0 0],nbcapteur,BC);
% Fart = Feta(:,:,iter_k) - T(k,tCalibration).*FDiskCentered;
% DiscreteDisk = DiscreteDiskCentered + tCalibration';
% [gStar,gTK,aTK] = HPenaltyOptimization(k,HEST(iter_k),Fart,Rin,tCalibration);

[gStar,gTK,aTK] = FArtPenalisation(Rin,tCalibration,BC,List_k(iter_k), Feta(:,:,iter_k), HEST(iter_k),DiscreteDiskCentered);
MHArt = HArtSoftDisk(k,Rin,tCalibration,nbcapteur);
H0gTK = MHArt*gTK;
H0gStar = MHArt*gStar;

%for iter_t = [3,11]
parfor iter_t = 1:Nt
    
    
    disp([num2str(iter_t) '/' num2str(Nt)]);
    t = [BX(iter_t),BY(iter_t)];
   
   [gStar,gTK,aTK] = FArtPenalisation(Rin,t,BC,List_k(iter_k), Feta(:,:,iter_k), HEST(iter_k),DiscreteDiskCentered);
   
    [MHArt,Nh] = HArtSoftDisk(k,Rin,t,nbcapteur);
    HgTK = MHArt*gTK;
    HgStar = MHArt*gStar;
    
    
    I1Temp(iter_t) = mean(vecnorm(imag(HgTK)));
    I2Temp(iter_t) = mean(vecnorm(imag(HgStar)));
    
    J1Temp(iter_t) = mean(vecnorm(HgTK - H0gTK));
    J2Temp(iter_t) = mean(vecnorm(HgStar - H0gStar));
    
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Inspection de la methode
     
%      load('RefPDESol/DnuRef.mat','DnuRef');
%      ThetaArt = 2*pi*(0:(1/Nh):(1-1/Nh));
%      SamplingPoint = 1;
%      z0 = DiscreteDiskCentered(SamplingPoint);
%      z = z0+t;
%      Hgzalpha = HgStar(:,SamplingPoint);
%        
%      
%      thetadisks = pi*(0:0.01:1.99);
%      ArtDisk = polyshape(Rin*cos(thetadisks)+BX(iter_t),Rin*sin(thetadisks)+BY(iter_t));
%      figure,
%      hold on 
%      plot(ArtDisk)
%      plot([Ax;Bx],[Ay;By],'red','LineWidth',3);
%      plot(z(1),z(2),'+black','LineWidth',2);
%      plot(BX(iter_t),BY(iter_t),'+black','LineWidth',2);
%      text(z(1)-.03,z(2)-.03,'z','FontSize',20)
%      text(BX(iter_t)+0.03,BY(iter_t)-0.03,'t_1','FontSize',20)
%      legend('$$\Omega_{t_1}$$','$$\Gamma$$','Interpreter','Latex','FontSize',25)
%      xlim([-.8 .7]);
%      ylim([-.8 .7]);
%      axis square
%      box('on');
%      if iter_t==3
%      print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/FixedSweepNoCrack.eps');
%      end 
%       if iter_t==11
%      print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/FixedSweepCrack.eps');
%      end  
%      
%      
%      linewidth = 3;
%      figure,hold on,
%      plot(ThetaArt,real(Hgzalpha),'LineWidth',linewidth);
%      plot(ThetaArt,real(DnuRef),'LineWidth',linewidth);
%      legend('$$\Re e \,H_{\partial \Omega}g_z^\alpha$$','$$\Re e \, \partial_\nu(w_z - \Phi_z)$$','Interpreter','Latex','FontSize',30)
%          xticks([0 pi/2 pi 3*pi/2 2*pi])   
%     xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
%     box('on');
%     if iter_t==3
%       print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RHgNoCrackR3k15.eps')
%     end
%     if iter_t==11
%       print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RHgCrackR3k15.eps')
%     end
%      
%       figure,hold on,
%      plot(ThetaArt,imag(Hgzalpha),'LineWidth',linewidth);
%      plot(ThetaArt,imag(DnuRef),'LineWidth',linewidth);
%      legend('$$\Im m \, H_{\partial \Omega}g_z^\alpha$$','$$\Im m \, \partial_\nu(w_z - \Phi_z)$$','Interpreter','Latex','FontSize',30)
%          xticks([0 pi/2 pi 3*pi/2 2*pi])   
%     xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
%     ylim([-2 2])
%     box('on');
%       if iter_t==3
%       print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/IHgNoCrackR3k15.eps')
%     end
%     if iter_t==11
%       print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/IHgCrackR3k15.eps')
%     end
     
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

    I1(:,:,iter_k) = I1Temp;
    I2(:,:,iter_k) = I2Temp;
    J1(:,:,iter_k) = J1Temp;
    J2(:,:,iter_k) = J2Temp;
    
    toc
    
    
   
 FileOutput = ['MultiFrqData/Monotonicity/FixedFrqMethod' num2str(iterNbCracks) '.mat'];
  save(FileOutput,'I1','J1','I2','J2','Rin','BXmin','BXmax','BYmin','BYmax','DB')

end 

end

%%
if DataSaveFlag==1
    save([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'I1','J1','I2','J2',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','-append')
end


 save([Folder 'FFRin01'  '.mat'],'I1','J1','I2','J2',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','-append')
    
    
%% Reload
tic
Folder = 'FixedFrqData/random2/';
 load([Folder 'FF.mat'],'FF','I1','J1','I2','J2',...
       'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By','List_k','I_LSM','I_GLSM')
   
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
Nk=length(List_k);

disp(['~~~~~~~~~~ ' Folder ' loaded' ' ~~~~~~~~~~']);
toc

 %% Regularization of the indicators

disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
RegList = sort([1:100 sqrt(2)*(1:100)]);

for RegLevel = 1
    RegLevel
    tic


% regularization of the previous indicators: see paper
I1Sum = zeros([size(BX) Nk]);I2Sum = zeros([size(BX) Nk]);J1Sum = zeros([size(BX) Nk]);
J2Sum = zeros([size(BX) Nk]);

% max of the indicators might be relevant
I1Max = zeros([size(BX) Nk]);I2Max = zeros([size(BX) Nk]);J1Max = zeros([size(BX) Nk]);
J2Max = zeros([size(BX) Nk]);

I1SumTemp = zeros(size(BX));
J1SumTemp = zeros(size(BX));
I2SumTemp = zeros(size(BX));
J2SumTemp = zeros(size(BX));
I1MaxTemp = zeros(size(BX));
I2MaxTemp = zeros(size(BX));
J1MaxTemp = zeros(size(BX));
J2MaxTemp = zeros(size(BX));



for iter_k = 1:Nk
    for iter_t = 1:Nt
    
    %disp([num2str(iter_t) '/' num2str(Nt)]);

    Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
    Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
    SBool = sum(sum(Bool));
       
    

        Temp = I1(:,:,iter_k);
        I1SumTemp(iter_t) = sum(Temp(Bool))/SBool;
        I1MaxTemp(iter_t) = max(Temp(Bool));

         Temp = I2(:,:,iter_k);
        I2SumTemp(iter_t) = sum(Temp(Bool))/SBool;
        I2MaxTemp(iter_t) = max(Temp(Bool));

        Temp = J1(:,:,iter_k);
        J1SumTemp(iter_t) = sum(Temp(Bool))/SBool;
        J1MaxTemp(iter_t) = max(Temp(Bool));

        Temp = J2(:,:,iter_k);
        J2SumTemp(iter_t) = sum(Temp(Bool))/SBool;
        J2MaxTemp(iter_t) = max(Temp(Bool));
        
  
    end
    
    I1Sum(:,:,iter_k) = I1SumTemp;
    I1Max(:,:,iter_k) = I1MaxTemp;
    
     I2Sum(:,:,iter_k) = I2SumTemp;
    I2Max(:,:,iter_k) = I2SumTemp;
    
    J1Sum(:,:,iter_k) = J1SumTemp;
    J1Max(:,:,iter_k) = J1MaxTemp;
    
    J2Sum(:,:,iter_k) = J2SumTemp;
    J2Max(:,:,iter_k) = J2MaxTemp;
end

toc
%

saveanothername([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],J2Sum,['J2SumReg' num2str(RegLevel)]);
saveanothername([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],I2Sum,['I2SumReg' num2str(RegLevel)]);
end

%%
save([Folder 'FF.mat'],'RegLevel','I1Sum','I2Sum','I1Max','I2Max','J1Sum','J2Sum','J1Max','J2Max','-append')
%% Plots


if Nk < 30
    
    for iter_k = 1:Nk
        
        
        FigIJ = figure('visible',PlotDisplayFlag);
        sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin)])

        subplot(2,2,1)
        hold on
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        %print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('I Tokhonov')

        subplot(2,2,2)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
       % print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('I GLSM')

        subplot(2,2,3)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J1(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('J Tikhonov')

        subplot(2,2,4)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('J GLSM')



        FigIJReg = figure('visible',PlotDisplayFlag);
        sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin) ', Reg = ' num2str(RegLevel)])

        subplot(2,4,1)
        hold on
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        %print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('ISum Tokhonov')

        subplot(2,4,2)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2Sum(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10,'MarkerFaceColor','white');
        axis square
        colormap(jet)
        colorbar
        %print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('ISum GLSM')

        subplot(2,4,3)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J1Sum(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('JSum Tikhonov')

        subplot(2,4,4)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('JSum GLSM')

        subplot(2,4,5)
        hold on
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Max(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        %print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('IMax Tokhonov')

        subplot(2,4,6)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2Max(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        %print(fig1,'-depsc2','/home/napal/nextcloud/WAVES2019/Images/temp1.eps')
        title('IMax GLSM')

        subplot(2,4,7)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J1Max(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('JMax Tikhonov')

        subplot(2,4,8)
        hold on, 
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Max(:,:,iter_k))
        set(gca,'Ydir','normal')
        plot([Ax Bx],[Ay By],'whiteX','MarkerSize',10);
        axis square
        colormap(jet)
        colorbar
        title('JMax GLSM')
    
    end 
end 

%% Save plots

if(PlotSaveFlag == 1)
    for iter_k = 1:Nk
        figtemp = figure, hold on,
        imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
        set(gca,'Ydir','normal') 
        plot([Ax Bx],[Ay By],'blacko','MarkerSize',7,'MarkerFaceColor','white');
        axis square
        colormap(jet)
        colorbar
        sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin)])
        print(figtemp,'-depsc2',[Folder num2str(List_k(iter_k)) 'k_JSum.eps'])
    end


    figControle = figure;%('visible','off');
    hold on, 
    %imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,rand(size((BX))))
    plot([Ax;Bx],[Ay;By],'red');
    set(gca,'Ydir','normal')
    axis square
    title('Geometry 5')
    print(figControle,'-depsc2',[Folder 'Geometry.eps'])
end

%% LSM - GLSM Indicators

ListX = BXmin:0.01:BXmax;
ListY = BYmin:0.01:BYmax;

I_LSM = zeros([length(ListX) length(ListY) Nk]);
I_GLSM = zeros([length(ListX) length(ListY) Nk]);

for iter_k = 3%:Nk

    disp([num2str(iter_k) '/' num2str(Nk)]);
    [I_LSM(:,:,iter_k),I_GLSM(:,:,iter_k)] = SamplingMethods(ListX,ListY,iter_k);
    
    
    
    

end

if(DataSaveFlag==1)
    save([Folder 'FF.mat'],'I_LSM','I_GLSM','ListX','ListY','-append');
end


%% LSM - GLSM ZOOM

BXminZoom = 0.5;
BXmaxZoom = 1.5;
BYminZoom = -1.5;
BYmaxZoom = -0.5;

figure, 
% cracks
hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);%,'Color', [.5 .5 .5]);
xlim([BXminZoom BXmaxZoom]);
ylim([BYminZoom BYmaxZoom]);
set(gca,'Ydir','normal')
axis square

ListXZoom = BXminZoom:0.001:BXmaxZoom;
ListYZoom = BYminZoom:0.001:BYmaxZoom;

I_LSMZoom = zeros([length(ListXZoom) length(ListYZoom) Nk]);
I_GLSMZoom = zeros([length(ListXZoom) length(ListYZoom) Nk]);

for iter_k = 3%:Nk

    disp([num2str(iter_k) '/' num2str(Nk)]);
    [I_LSMZoom(:,:,iter_k),I_GLSMZoom(:,:,iter_k)] = SamplingMethods(ListXZoom,ListYZoom,iter_k);

end

 save([Folder 'FF.mat'],'I_LSMZoom','I_GLSMZoom','ListXZoom','ListYZoom','-append');

%%

for iter_k = 1:Nk
figure, 
sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin)])
subplot(2,3,1)
hold on,
imagesc(ListX,ListY,I_LSM(:,:,iter_k).^(-1))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal') 
%plot([Ax Bx],[Ay By],'whiteo','MarkerSize',10);
axis square
colormap(jet)
colorbar

subplot(2,3,3)
hold on,
imagesc(ListX,ListY,I_GLSM(:,:,iter_k).^(-1))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal') 
%plot([Ax Bx],[Ay By],'whiteo','MarkerSize',10);
axis square
colormap(jet)
colorbar

subplot(2,3,2)
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
%plot([Ax Bx],[Ay By],'whiteo','MarkerSize',10);
axis square
colormap(jet)
colorbar
title('JSum GLSM')

subplot(2,3,5)
hold on, 
plot([Ax;Bx],[Ay;By],'red');
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
axis square
colormap('jet');
colorbar;
title('Resolution of the image');

end 






        
     
