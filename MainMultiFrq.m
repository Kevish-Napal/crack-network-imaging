set(0,'DefaultFigureWindowStyle','docked')

% Cleaning
clear all
close all
clc


addpath('FGen')

%% Parameters

%global DiscreteDiskCentered Feta HEST List_k % discretization of the centered artificial disk

% Choose Geometry between the available geometries:
%Folder = 'MultiFrqData/1Crack/';%'Geometry_6Cracks_k14.mat';%'Geometry_2Cracks_k14.mat';


Rin = 0.1; 

for iterNbCracks = 3:9
    iterNbCracks
    
FileOutput = ['MultiFrqData/Monotonicity/FF' num2str(iterNbCracks) '.mat'];
load(FileOutput,'FF','Ax','Ay','Bx','By','List_k');
%load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'Ax','Ay','Bx','By','FF','List_k');

[nbcapteur,~,Nk] = size(FF); 


% redimensionnement de List_k
% List_k = List_k(1:10:Nk)
% FF = FF(:,:,1:10:Nk);
% Nk = length(List_k);


% Choose resolution of the scan


BXmin = -1.; % range of the centers coordinates
BXmax = 1.;
BYmin = -1.;
BYmax = 1.;
DB = Rin; % distance between each centers
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
dk = List_k(2)-List_k(1);

% Choose artificial background boundary condition:
% BC1 u + BC2 d_n u = 0.
BC = [1 0];

% plot and/or save computed indicators?
PlotDisplayFlag = 0; 
DataSaveFlag = 1;
PlotSaveFlag = 0;

% Choose noise level of farfield in the next section (eta).

% Control of the scan parameters
figControle = figure;%('visible','off');
clearvars('figProbe','figCracks')
hold on, 
figCracks = plot([Ax;Bx],[Ay;By],'red','LineWidth',5);
figProbe = plot([BXmin;BXmax;BXmax;BXmin;BXmin],[BYmin;BYmin;BYmax;BYmax;BYmin],'blue','LineWidth',5);
set(gca,'Ydir','normal')
axis square
title('Resolution of the image');
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',15)
xlim(1.2*[BXmin BXmax]);
ylim(1.2*[BYmin BYmax]);


% visualisation of the artificial backgrounds (when possible):
if Nt<1500
poly = DiskSwipe(BX,BY,Rin,30);
figure,hold on,
plot(poly,'FaceAlpha',0.1,'EdgeAlpha',0);
plot([Ax;Bx],[Ay;By],'red','LineWidth',1);
axis equal
title('The artificial backgrounds');
end

% dirichlet eigenvalue
kmax = 10./Rin;
k = 0:0.001:kmax;

tic
[DE,pmax] = DirichletEigenvalues(kmax*Rin);
toc

figure,hold on,
for p = 0:pmax

    BJpkR = besselj(p,k*Rin);
    figBessel = plot(k,BJpkR,'--');
    figBessel.Color = [220,220,220]*0.001;
end


xlabel('k');
ylabel('J_p(k*Rin)')

DE = 1/Rin*DE;
figDE = plot(DE,zeros(size(DE)),'ro');
figList_k = plot(List_k,zeros(size(List_k)),'blue+');
legend([figBessel figDE figList_k],'J_p(k*Rin)','Dirichlet Eigenvalues','Sampling k');
title(['Rin = ' num2str(Rin)]);

%% Add of noise on the data

eta=0.015;
Feta = zeros(nbcapteur,nbcapteur,Nk);
HEST = zeros(1,Nk);
nF = zeros(1,Nk);

for i=1:Nk

     
    % Bruit uniforme
    % noise=1+eta*(rand(size(FF))*2-1)+sqrt(-1)*eta*(rand(size(FF))*2-1);

    % Bruit Gaussien
    noise=1+eta*(randn(size(FF(:,:,i))))+sqrt(-1)*eta*(randn(size(FF(:,:,i))));

    Feta(:,:,i)=FF(:,:,i).*noise;

    HEST(i)=norm(Feta(:,:,i)-FF(:,:,i));
    if HEST(i) < 1.e-02
        HEST(i)= 1.e-02;                       
    end
    
    nF(i)=norm(Feta(:,:,i));
    %disp(['Bruit en pourcentage = ',num2str(HEST(i)/nF(i)*100)])
end

NoisePercentageLevel = mean(100*HEST./nF)

%% Discretization of the disk

% Sampling points on the centered artificial background
Nbz = 300; % Nb Approx de points de discretisation interieurs du disque artificiel
[X,Y]=meshgrid(linspace(-Rin,Rin,sqrt(4./pi*Nbz)));
Bool = X.^2+Y.^2 < (Rin).^2;
DiscreteDiskCentered = transpose([X(Bool),Y(Bool)]);
NbPoints = length(DiscreteDiskCentered); % Nb de points de discretisation interieurs du disque artificiel
figure,plot(X(Bool),Y(Bool),'red+'), 
hold on, plot(Rin*cos(0:0.1:7),Rin*sin(0:0.1:7),'black'),
axis square;
title(['Discretisation interieur du disque artificiel centre de rayon ' num2str(Rin) ', nbPtsInterieur: ' num2str(NbPoints)]);
                                                                                                             
%%

Scalaire1Temp = zeros(1,Nt);
Scalaire2Temp = zeros(1,Nt);
Scalaire3Temp = zeros(1,Nt);
Scalaire = zeros(Nk,Nt,3);
        
for iter_k = 1:Nk % Boucles sur k   
    tic
    disp(['k= ' num2str(iter_k) '/' num2str(Nk)]);
    k = List_k(iter_k);
    h = HEST(iter_k);
    FTemp = Feta(:,:,iter_k);
    parfor iter_t = 1:Nt

     % disp(['t= ' num2str(iter_t) '/' num2str(Nt) ' ,k= ' num2str(iter_k) '/' num2str(Nk)]);
        
        t = [BX(iter_t),BY(iter_t)];

       [gStar,gTK,aTK,Penalisation] = FArtPenalisation(Rin,t,BC,k,FTemp,h,DiscreteDiskCentered);

        Scalaire1Temp(iter_t) = norm(norm(HArtSoftDisk(k,Rin,t,nbcapteur)*gTK));
        Scalaire2Temp(iter_t) = norm(norm(HArtSoftDisk(k,Rin,t,nbcapteur)*gStar));
        Scalaire3Temp(iter_t) = norm(dot(Penalisation*gStar,gStar));
         
         
    end

    Scalaire(iter_k,:,1) = Scalaire1Temp;
    Scalaire(iter_k,:,2) = Scalaire2Temp;
    Scalaire(iter_k,:,3) = Scalaire3Temp;

    toc
end
    

if DataSaveFlag==1

save(FileOutput,'Scalaire',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','-append')
end

   
end
%% load existing data fo plots
% Folder = 'FullGeometry/';
%  load([Folder 'FF.mat'],'I1','J1','I2','J2',...
%        'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By','List_k')
%    
% [BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
% Nt = numel(BX);
% Nk=length(List_k);

%% Controle des resultats
Rin = 0.5
Folder = 'MultiFrqData/random/';
load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'Ax','Ay','Bx','By','FF','List_k',...
'Scalaire','Rin','BXmin','BXmax','BYmin','BYmax','DB');

Nt = numel(BX);

%%
tic

DE = 1/Rin*DirichletEigenvalues(List_k(end)*Rin);
Density1 = zeros(size(BX));
Density2 = zeros(size(BX));
Density3 = zeros(size(BX));
PKS = {};
ScalaireGolay = zeros(size(Scalaire));

for iter_t = 1:Nt 
    
    ScalaireGolay(:,iter_t,1) = sgolayfilt(Scalaire(:,iter_t,1),1,3);
    ScalaireGolay(:,iter_t,2) = sgolayfilt(Scalaire(:,iter_t,2),1,3);
    ScalaireGolay(:,iter_t,3) = sgolayfilt(Scalaire(:,iter_t,3),1,3);

    [pks,I] = findpeaks(ScalaireGolay(:,iter_t,1),'MinPeakProminence',3*dk);
    Density1(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{1,iter_t} = [pks,I];
    
    
    [pks,I] = findpeaks(ScalaireGolay(:,iter_t,2),'MinPeakProminence',3*dk);
    Density2(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{2,iter_t} = [pks,I];
    
    [pks,I] = findpeaks(ScalaireGolay(:,iter_t,3),'MinPeakProminence',3*dk);
    Density3(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{3,iter_t} = [pks,I];
    

end

toc
%%



% visualisation of the artificial backgrounds (when possible):
if Nt<1500
poly = DiskSwipe(BX,BY,Rin,30);
figure,hold on,
plot(poly,'FaceAlpha',0.1,'EdgeAlpha',0);
plot([Ax;Bx],[Ay;By],'red','LineWidth',1);
axis equal
title('The artificial backgrounds');
end

fig1 = figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density1))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'red');
axis equal
colormap(jet)
colorbar



fig2 = figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density2))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'red');
axis equal
colorbar
colormap(jet)



fig3 = figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density3))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'red');
axis equal
colorbar
colormap(jet)

thetadisks = pi*(0:0.01:1.99);
% les 5 indicatrices pour le titre des affichages
S = {'$$\langle \widetilde{F}_\# g_\#,g_\#\rangle$$',...
'$$\|\widetilde{H}g_{TK}\|$$',...
'$$\|H_{art}g_{TK}\|$$',...
'$$\|\widetilde{H}g_\#\|$$',...
'$$\|H_{art}g_\#\|$$'};

%%
Klim = 10.;
K = 0:0.01:Klim;
bool = List_k < Klim; 
itersave = 0;
 indic = 3;
for iter_t = 1:10%numel(BX)
 
    ArtDisk = polyshape(Rin*cos(thetadisks)+BX(iter_t),Rin*sin(thetadisks)+BY(iter_t));
    figure,
    pos1 = [-0.02 0.3 0.3 0.3];
    pos2 = [0.28 0.1 0.7 0.8];
    subplot('Position',pos2)
    hold on,
    plot(List_k(:),ScalaireGolay(:,iter_t,indic),'--blue')
    
    pks = PKS{indic,iter_t}(:,1);
    I = PKS{indic,iter_t}(:,2);
    plot(List_k(I),pks,'r+')
    for de = 1:2
       xline(DE(de));
    end

    subplot('Position',pos1)
    hold on 
    plot(ArtDisk)
    plot([Ax;Bx],[Ay;By],'red','LineWidth',1);
    xlim([BXmin BXmax]);
    ylim([BYmin BYmax]);
    axis square
    title(['t = ' num2str([BX(iter_t),BY(iter_t)]) ', ' S{indic}],'interpreter','latex')
    legend('g')

  
end


 %%
    
disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
tic
RegLevel = 1
RegList = sort([1:10 sqrt(2)*(1:10)]);

% regularization of the previous indicators: see paper
I1Sum = zeros([size(BX)]);

% max of the indicators might be relevant
I1Max = zeros([size(BX)]);

I1SumTemp = zeros(size(BX));




for iter_k = 1:Nk
    for iter_t = 1:Nt
    
    disp([num2str(iter_t) '/' num2str(Nt)]);

    Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
    Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
    SBool = sum(sum(Bool));
       
        
        I1Sum(iter_t) = sum(Density3(Bool))/SBool;
        I1Max(iter_t) = max(Density3(Bool));

        
  
    end
end

toc

%%

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'red');
axis equal
colormap(jet)
colorbar

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Max)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'red');
axis equal
colormap(jet)
colorbar
