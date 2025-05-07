set(0,'DefaultFigureWindowStyle','docked')
close all
clear all
clc

PlotSaveFlag = 0;


Color.white = [1 1 1];
Color.black = [0 0 0];
Color.gray1 = [0.85 0.85 0.85];
Color.gray2 = [0.7 0.7 0.7];
Color.gray3 = [0.4 0.4 0.4];
Color.gray4 = [0.2 0.2 0.2];
Color.mystic = [219, 76, 119]/256;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Multi Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% load Random cracks
Folder = 'MultiFrqData/random/';%'Geometry_6Cracks_k14.mat';%'Geometry_2Cracks_k14.mat';
Rin = 0.1; 
load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'Ax','Ay','Bx','By','List_k',...
'Scalaire','Rin','BXmin','BXmax','BYmin','BYmax','DB');


% Plot of the crack distribution

Gauss = @(x,mu,sigma) 1/(sigma*sqrt(2*pi))*exp(-0.5*((x-mu)/sigma).^2)

dx = -2.5:0.001:2.5;
[X,Y] = meshgrid(dx);

s1 = 1;
s2 = 0.5;
s3 = .4;
s4 = .2;

CrackDensity =  Gauss(X,-1,s1).*Gauss(Y,1,s1) ...
              + Gauss(X,1,s2).*Gauss(Y,1,s2) ...;
              + Gauss(X,-1,s3).*Gauss(Y,-1,s3) ...; 
              + Gauss(X,1,s4).*Gauss(Y,-1,s4);  

          
%% Computation of indicators from Scalaire
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
Nk = length(List_k);
dk = List_k(2)-List_k(1)
tic

DE = 1/Rin*DirichletEigenvalues(List_k(end)*Rin);
Density1 = zeros(size(BX));
Density2 = zeros(size(BX));
Density3 = zeros(size(BX));
PKS = {};

for iter_t = 1:Nt 
    
    Scalaire(:,iter_t,1) = sgolayfilt(Scalaire(:,iter_t,1),1,15);
    Scalaire(:,iter_t,2) = sgolayfilt(Scalaire(:,iter_t,2),1,15);
    Scalaire(:,iter_t,3) = sgolayfilt(Scalaire(:,iter_t,3),1,15);

    [pks,I] = findpeaks(Scalaire(:,iter_t,1),'MinPeakProminence',3*dk);
    Density1(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{1,iter_t} = [pks,I];
    
    
    [pks,I] = findpeaks(Scalaire(:,iter_t,2),'MinPeakProminence',3*dk);
    Density2(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{2,iter_t} = [pks,I];
    
    [pks,I] = findpeaks(Scalaire(:,iter_t,3),'MinPeakProminence',3*dk);
    Density3(iter_t) = sum(min(abs(List_k(I) - DE')));
    PKS{3,iter_t} = [pks,I];
    

end

toc

%% Regularizations
disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
tic
RegLevel = 5
RegList = sort([1:10 sqrt(2)*(1:10)]);

% regularization of the previous indicators: see paper
I1Sum = zeros([size(BX)]);

% max of the indicators might be relevant
I1Max = zeros([size(BX)]);

I1SumTemp = zeros(size(BX));




for iter_k = 1:Nk
    for iter_t = 1:Nt
    
    %disp([num2str(iter_t) '/' num2str(Nt)]);

    Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
    Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
    SBool = sum(sum(Bool));
       
        
        I1Sum(iter_t) = sum(Density3(Bool))/SBool;
        I1Max(iter_t) = max(Density3(Bool));

        
  
    end
end

toc

% Plots
Crackeps = 3;
CrackColor = Color.gray2;

theta = 0.01:0.01:2*pi;


    
figure,hold on 
%imagesc(dx,dx,CrackDensity/4)
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', Color.black);
figProbe = plot([BXmin;BXmax;BXmax;BXmin;BXmin],[BYmin;BYmin;BYmax;BYmax;BYmin],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',25)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/CrackDistribution.eps')

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color',CrackColor);
axis equal
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomMultR1.eps')

figure;
hold on, 
contourf(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color',CrackColor);
axis equal
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
colormap(jet)
colorbar


%% load Straight cracks
tic
load('../TestFF2019/CrackDensities.mat','Density3')
toc

Rin = 1;
%Zone de balayage:
BXmin = 1;
BXmax = 26;
BYmin = 1;
BYmax = 26;
DB = 1; % pas du balayage
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX)

% Les fissures
CentreX = [-2,  2-0.3, 2+0.3,    -2-0.2, -2, -2+0.2,    2-0.4, 2-0.2, 2, 2+0.2, 2+0.4];
CentreY = [ 2,  2    , 2    ,    -2    , -2, -2    ,    -2   ,-2    ,-2,  -2  ,   -2];
r = 0.25;
alpha = pi/2.;
cosinus = r/2*cos(alpha);
sinus = r/2*sin(alpha);
Ax = cosinus + CentreX;  
Ay = sinus + CentreY;
Bx = -cosinus + CentreX;
By = -sinus + CentreY;


figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,Density3)
set(gca,'Ydir','normal')
axis equal
colormap(jet)
colorbar


%% Regularizations
disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
tic
RegList = sort([1:10 sqrt(2)*(1:10)]);
for RegLevel = 2


% regularization of the previous indicators: see paper
I1Sum = zeros([size(BX)]);

% max of the indicators might be relevant
I1Max = zeros([size(BX)]);

I1SumTemp = zeros(size(BX));




for iter_k = 1:Nk
    for iter_t = 1:Nt
    
    %disp([num2str(iter_t) '/' num2str(Nt)]);

    Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
    Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
    SBool = sum(sum(Bool));
       
        
        I1Sum(iter_t) = sum(Density3(Bool))/SBool;
        I1Max(iter_t) = max(Density3(Bool));

        
  
    end
end

toc

% Plots
Crackeps = 3;

theta = 0.01:0.01:2*pi;

Nmap = 1000;
mymap = zeros(Nmap,3);
mymap(1,:) = 1;
mymap(2:end,1) = 1;
mymap(2:end,2) = 1-(2:Nmap)/Nmap;
mymap(2:end,3) = 1-(2:Nmap)/Nmap;

figure;
hold on, 
imagesc(-3:3,-3:3,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
axis equal
colormap(jet)
colorbar
xlim([-3 3]);
ylim([-3 3]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/StraightMultR5Reg2.eps')


end


%% Single Crack
Folder = 'MultiFrqData/1Crack/';%'Geometry_6Cracks_k14.mat';%'Geometry_2Cracks_k14.mat';
Rin = 0.3; 
load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'Ax','Ay','Bx','By','List_k',...
'Scalaire','Rin','BXmin','BXmax','BYmin','BYmax','DB');


Crackeps = 3;
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
dk = List_k(2) - List_k(1);
% visualisation of the artificial backgrounds (when possible):
if Nt<1500
poly = DiskSwipe(BX,BY,Rin,30);
figure,hold on,
plot([Ax;Bx],[Ay;By],'red','LineWidth',3);
plot(poly,'FaceAlpha',0.1,'EdgeAlpha',0);
legend('$$\Gamma$$','$$\Omega_t$$','Interpreter','Latex','FontSize',25)
axis equal
xlim([-.8 .7]);
ylim([-.8 .7]);
end
box('on')
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultSweep.eps')

tic

DE = 1/Rin*DirichletEigenvalues(List_k(end)*Rin);
Density1 = zeros(size(BX));
Density2 = zeros(size(BX));
Density3 = zeros(size(BX));
PKS = {};
ScalaireGolay = zeros(size(Scalaire));

thetadisks = pi*(0:0.01:1.99);

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

indic = 3;
iter_t = 3

ArtDisk = polyshape(Rin*cos(thetadisks)+BX(iter_t),Rin*sin(thetadisks)+BY(iter_t));


figure,
hold on,
plot(List_k(:),ScalaireGolay(:,iter_t,indic),'--blue')

pks = PKS{indic,iter_t}(:,1);
I = PKS{indic,iter_t}(:,2);
plot(List_k(I),pks,'ro','MarkerSize',7,'MarkerFaceColor','red')
for de = DE(1)
   xline(de);
end
xlim([1 10]);
box('on')
xlabel('$$k$$','Interpreter','Latex','FontSize',30)
ylabel('$$\int_{\Omega_{t_1}} P(g_z^n(k))\,dz$$','Interpreter','Latex','FontSize',30,'Position',[.9 8])

print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/ScalaireNoCrack.eps')

figure,
hold on 
plot(ArtDisk)
plot([Ax;Bx],[Ay;By],'red','LineWidth',3);
plot(BX(iter_t),BY(iter_t),'black+','LineWidth',3)
text(BX(iter_t),BY(iter_t)-.04,'t_1','FontSize',20)
xlim([-.8 .7]);
ylim([-.8 .7]);
legend('$$\Omega_{t_1}$$','$$\Gamma$$','Interpreter','Latex','FontSize',25)
axis square
box('on')
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/MultSweepNoCrack.eps')

iter_t = 11

ArtDisk = polyshape(Rin*cos(thetadisks)+BX(iter_t),Rin*sin(thetadisks)+BY(iter_t));
figure,
hold on,
plot(List_k(:),ScalaireGolay(:,iter_t,indic),'--blue')

pks = PKS{indic,iter_t}(:,1);
I = PKS{indic,iter_t}(:,2);
plot(List_k(I),pks,'ro','MarkerSize',7,'MarkerFaceColor','red')
for de = DE(1)
xline(de);
end
xlim([1 10]);
box('on')
xlabel('$$k$$','Interpreter','Latex','FontSize',30)
ylabel('$$\int_{\Omega_{t_2}} P(g_z^n(k))\,dz$$','Interpreter','Latex','FontSize',30,'Position',[.9 80])
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/ScalaireCrack.eps')

figure,
hold on 
plot(ArtDisk)
plot([Ax;Bx],[Ay;By],'red','LineWidth',3);
plot(BX(iter_t),BY(iter_t),'black+','LineWidth',3)
text(BX(iter_t),BY(iter_t)-.04,'t_2','FontSize',20)
legend('$$\Omega_{t_2}$$','$$\Gamma$$','Interpreter','Latex','FontSize',25)
xlim([-.8 .7]);
ylim([-.8 .7]);
axis square
box('on')
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/MultSweepCrack.eps')


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


CrackColor = Color.gray2;

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density3))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis equal
colorbar
colormap(jet)
xlim([-0.65 0.55]);
ylim([-0.65 0.55]);

%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultR3.eps')

toc
  
%% Single crack DEmo HR
Folder = 'MultiFrqData/1Crack/';%'Geometry_6Cracks_k14.mat';%'Geometry_2Cracks_k14.mat';
Rin = 0.1; 
load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'Ax','Ay','Bx','By','List_k',...
'Scalaire','Rin','BXmin','BXmax','BYmin','BYmax','DB');


Crackeps = 3;
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
dk = List_k(2) - List_k(1);


tic

DE = 1/Rin*DirichletEigenvalues(List_k(end)*Rin);
Density1 = zeros(size(BX));
Density2 = zeros(size(BX));
Density3 = zeros(size(BX));
PKS = {};
ScalaireGolay = zeros(size(Scalaire));

thetadisks = pi*(0:0.01:1.99);

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

indic = 3;


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


CrackColor = Color.gray2;

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density3))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis equal
colorbar
colormap(jet)
xlim([-0.65 0.55]);
ylim([-0.65 0.55]);

print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultR1.eps')
toc

%% Regularization
disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
tic
RegList = sort([1:10 sqrt(2)*(1:10)]);

for RegLevel = 1



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

figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis equal
colorbar
colormap(jet)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);

print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultR1Reg1.eps')


toc

end

%% Monotonicity %%%%

global Feta HEST List_k

for iterNbCracks = 1:8

FileOutput = ['MultiFrqData/Monotonicity/FF' num2str(iterNbCracks) '.mat'];
load(FileOutput,'Ax','Ay','Bx','By','List_k',...
'Scalaire','Rin','BXmin','BXmax','BYmin','BYmax','DB','FF');


% Cracks distribution

% Fissures distribution zoom

figure,hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
figProbe = plot([BXmin;BXmax;BXmax;BXmin;BXmin],[BYmin;BYmin;BYmax;BYmax;BYmin],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',25,'Location','northwest')
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square
print(gcf,'-depsc2',['/Users/kev/Downloads/version7/Images/MonotonyDistrib'  num2str(iterNbCracks) '.eps'])





Crackeps = 4;
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
dk = List_k(2) - List_k(1);

% % visualisation of the artificial backgrounds (when possible):
% if Nt<1500
% poly = DiskSwipe(BX,BY,Rin,30);
% figure,hold on,
% plot([Ax;Bx],[Ay;By],'red','LineWidth',3);
% plot(poly,'FaceAlpha',0.1,'EdgeAlpha',0);
% legend('$$\Gamma$$','$$\Omega_t$$','Interpreter','Latex','FontSize',25)
% axis equal
% xlim([BXmin BXmax]);
% ylim([BYmin BYmax]);
% end
% box('on')
% % print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultSweep.eps')

tic

DE = 1/Rin*DirichletEigenvalues(List_k(end)*Rin);
Density1 = zeros(size(BX));
Density2 = zeros(size(BX));
Density3 = zeros(size(BX));
PKS = {};
ScalaireGolay = zeros(size(Scalaire));

thetadisks = pi*(0:0.01:1.99);

for iter_t = 1:Nt 
    
    ScalaireGolay(:,iter_t,1) = sgolayfilt(Scalaire(:,iter_t,1),1,11);
    ScalaireGolay(:,iter_t,2) = sgolayfilt(Scalaire(:,iter_t,2),1,11);
    ScalaireGolay(:,iter_t,3) = sgolayfilt(Scalaire(:,iter_t,3),1,11);

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

%% Control
% thetadisks = pi*(0:0.01:1.99);
% % les 5 indicatrices pour le titre des affichages
% S = {'$$\langle \widetilde{F}_\# g_\#,g_\#\rangle$$',...
% '$$\|\widetilde{H}g_{TK}\|$$',...
% '$$\|H_{art}g_{TK}\|$$',...
% '$$\|\widetilde{H}g_\#\|$$',...
% '$$\|H_{art}g_\#\|$$'};
% 
% 
% Klim = 10.;
% K = 0:0.01:Klim;
% bool = List_k < Klim; 
% itersave = 0;
%  indic = 3;
% for iter_t = 220:270%numel(BX)
%  
%     ArtDisk = polyshape(Rin*cos(thetadisks)+BX(iter_t),Rin*sin(thetadisks)+BY(iter_t));
%     figure,
%     pos1 = [-0.02 0.3 0.3 0.3];
%     pos2 = [0.28 0.1 0.7 0.8];
%     subplot('Position',pos2)
%     hold on,
%     plot(List_k(:),ScalaireGolay(:,iter_t,indic),'--blue')
%     
%     pks = PKS{indic,iter_t}(:,1);
%     I = PKS{indic,iter_t}(:,2);
%     plot(List_k(I),pks,'r+')
%     for de = 1:2
%        xline(DE(de));
%     end
% 
%     subplot('Position',pos1)
%     hold on 
%     plot(ArtDisk)
%     plot([Ax;Bx],[Ay;By],'red','LineWidth',1);
%     xlim([BXmin BXmax]);
%     ylim([BYmin BYmax]);
%     axis square
%     title(['t = ' num2str([BX(iter_t),BY(iter_t)]) ', ' S{indic}],'interpreter','latex')
%     legend('g') 
% end

%%



CrackColor = Color.gray2;

% figure;
% hold on, 
% imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,abs(Density3))
% set(gca,'Ydir','normal')
% plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
% axis equal
% colorbar
% colormap(jet)
% xlim([BXmin BXmax]);
% ylim([BYmin BYmax]);
% %print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoMultR3.eps')

% Regularization
disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
tic
RegList = sort([1:10 sqrt(2)*(1:10)]);
RegLevel = 1
% regularization of the previous indicators: see paper
I1Sum = zeros([size(BX)]);
% max of the indicators might be relevant
I1Max = zeros([size(BX)]);

    for iter_t = 1:Nt  
    disp([num2str(iter_t) '/' num2str(Nt)]);
    Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
    Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
    SBool = sum(sum(Bool));        
    I1Sum(iter_t) = sum(Density3(Bool))/SBool;
    I1Max(iter_t) = max(Density3(Bool));
    end


figure;
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I1Sum)
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis equal
colorbar
colormap(jet)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2',['/Users/kev/Downloads/version7/Images/MonotonyMult'  num2str(iterNbCracks) '.eps'])


disp(['~~~~~~~~~~ FACTORIZATION METHOD ~~~~~~~~~~']);

[nbcapteur,~,Nk] = size(FF);

eta=.02;
Feta = zeros(nbcapteur,nbcapteur,Nk);
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

% Factorisations method


iter_k = 251
    ListX = BXmin:.01:BXmax;
    ListY = BYmin:.01:BYmax;
 Mreg = 100

I_FM = CrackFM(ListX,ListY,iter_k,Mreg);
figure,
hold on,
imagesc(ListX,ListY,I_FM)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2',['/Users/kev/Downloads/version7/Images/MonotonyFM'  num2str(iterNbCracks) '.eps'])



% load for fixed frq method

FileOutput = ['MultiFrqData/Monotonicity/FixedFrqMethod' num2str(iterNbCracks) '.mat'];
 load(FileOutput,'I1','J1','I2','J2','Rin','BXmin','BXmax','BYmin','BYmax','DB');
   
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
Nk=length(List_k);

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
set(gca,'Ydir','normal')
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2',['/Users/kev/Downloads/version7/Images/MonotonyFixed'  num2str(iterNbCracks) '.eps'])
        
        
% Regularization of the indicators
% 
% disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
% RegList = sort([1:100 sqrt(2)*(1:100)]);

% for RegLevel = 1
%     RegLevel
%     
% 
% 
% % regularization of the previous indicators: see paper
% I1Sum = zeros([size(BX) Nk]);I2Sum = zeros([size(BX) Nk]);J1Sum = zeros([size(BX) Nk]);
% J2Sum = zeros([size(BX) Nk]);
% 
% 
% I1SumTemp = zeros(size(BX));
% J1SumTemp = zeros(size(BX));
% I2SumTemp = zeros(size(BX));
% J2SumTemp = zeros(size(BX));



% for iter_k = Nk
%     for iter_t = 1:Nt
%     
%     disp([num2str(iter_t) '/' num2str(Nt)]);
% 
%     Alentours = sqrt((BX - BX(iter_t)).^2 + (BY - BY(iter_t)).^2);
%     Bool = Alentours < (RegList(RegLevel)+0.2)*Rin;
%     SBool = sum(sum(Bool));
%        
%     
% 
%         Temp = I1(:,:,iter_k);
%         I1SumTemp(iter_t) = sum(Temp(Bool))/SBool;
%         I1MaxTemp(iter_t) = max(Temp(Bool));
% 
%          Temp = I2(:,:,iter_k);
%         I2SumTemp(iter_t) = sum(Temp(Bool))/SBool;
%         I2MaxTemp(iter_t) = max(Temp(Bool));
% 
%         Temp = J1(:,:,iter_k);
%         J1SumTemp(iter_t) = sum(Temp(Bool))/SBool;
%         J1MaxTemp(iter_t) = max(Temp(Bool));
% 
%         Temp = J2(:,:,iter_k);
%         J2SumTemp(iter_t) = sum(Temp(Bool))/SBool;
%         J2MaxTemp(iter_t) = max(Temp(Bool));
%         
%   
%     end
%     
%     I1Sum(:,:,iter_k) = I1SumTemp;
%     I1Max(:,:,iter_k) = I1MaxTemp;
%     
%      I2Sum(:,:,iter_k) = I2SumTemp;
%     I2Max(:,:,iter_k) = I2SumTemp;
%     
%     J1Sum(:,:,iter_k) = J1SumTemp;
%     J1Max(:,:,iter_k) = J1MaxTemp;
%     
%     J2Sum(:,:,iter_k) = J2SumTemp;
%     J2Max(:,:,iter_k) = J2MaxTemp;
% end
% 
% toc
% end
% 
% FigIJ = figure('visible',PlotDisplayFlag);
%         sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin)])
% 
% 
% 
%         subplot(1,2,1)
%         hold on, 
%         imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J1(:,:,iter_k))
%         set(gca,'Ydir','normal')
% %         plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
%         axis square
%         colormap(jet)
%         colorbar
%         title('J Tikhonov')
%         xlim([BXmin BXmax]);
%         ylim([BYmin BYmax]);
% 
%         subplot(1,2,2)
%         hold on, 
%         imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
%         set(gca,'Ydir','normal')
% %         plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
%         axis square
%         colormap(jet)
%         colorbar
%         title('J GLSM')
%         xlim([BXmin BXmax]);
%         ylim([BYmin BYmax]);
% 
% 
% 
%         FigIJReg = figure('visible',PlotDisplayFlag);
%         sgtitle(['lambda = ' num2str(2*pi/List_k(iter_k)) ', k = ' num2str(List_k(iter_k)) ', r = ' num2str(Rin) ', Reg = ' num2str(RegLevel)])
% 
%         subplot(1,2,1)
%         hold on, 
%         imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J1Sum(:,:,iter_k))
%         set(gca,'Ydir','normal')
% %         plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
%         axis square
%         colormap(jet)
%         colorbar
%         title('JSum Tikhonov')
%         xlim([BXmin BXmax]);
%         ylim([BYmin BYmax]);
% 
%         subplot(1,2,2)
%         hold on, 
%         imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
%         set(gca,'Ydir','normal')
% %         plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
%         axis square
%         colormap(jet)
%         colorbar
%         title('JSum GLSM')
%         xlim([BXmin BXmax]);
%         ylim([BYmin BYmax]);
%     

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%% Rin = 0.01
tic
Folder = 'FixedFrqData/random2/';
 load([Folder 'FF.mat'],'FF','I1','J1','I2','J2',...
       'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By','List_k','I_LSM','I_GLSM','ListX','ListY',...
       'RegLevel','I1Sum','I2Sum','I1Max','I2Max','J1Sum','J2Sum','J1Max','J2Max',...
       'I2SumReg25','I2SumReg35','I2SumReg40','I2SumReg50',...
       'J2SumReg25','J2SumReg35','J2SumReg40','J2SumReg50',...
       'I_LSMZoom','I_GLSMZoom','ListXZoom','ListYZoom')
   
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nt = numel(BX);
Nk=length(List_k);
disp(['~~~~~~~~~~ ' Folder ' loaded' ' ~~~~~~~~~~']);
disp(['Rin = ' num2str(Rin)]);
disp(['RegLevel = ' num2str(RegLevel)]);
toc

%%


iter_k = 3

CrackColor = Color.gray2;

Crackeps = 3;

figure, 
% cracks
hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);%,'Color', [.5 .5 .5]);
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square

% GLSM
figure,
hold on,
imagesc(ListX,ListY,I_GLSM(:,:,iter_k).^(-1))
xlim([-2.5 2.5]);
ylim([-2 2.5]);
set(gca,'Ydir','normal') 
%plot([Ax Bx],[Ay By],'whiteo','MarkerSize',10);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomGLSMk20.eps')

% Fissures distribution zoom
BXminZoom = ListXZoom(1);
BXmaxZoom =  ListXZoom(end);
BYminZoom = ListYZoom(1);
BYmaxZoom =  ListYZoom(end);
figure,hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
figProbe = plot([BXminZoom;BXmaxZoom;BXmaxZoom;BXminZoom;BXminZoom],[BYminZoom;BYminZoom;BYmaxZoom;BYmaxZoom;BYminZoom],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',25)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/CrackDistributionZoom.eps')

figure,hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
figProbe = plot([BXminZoom;BXmaxZoom;BXmaxZoom;BXminZoom;BXminZoom],[BYminZoom;BYminZoom;BYmaxZoom;BYmaxZoom;BYminZoom],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',15)
xlim([BXminZoom BXmaxZoom]);
ylim([BYminZoom BYmaxZoom]);
set(gca,'Ydir','normal')
axis square
grid on
grid minor



% GLSMZoom
figure,
hold on,
imagesc(ListXZoom,ListYZoom,I_GLSMZoom(:,:,iter_k).^(-1))
xlim([ListXZoom(1) ListXZoom(end)]);
ylim([ListYZoom(1) ListYZoom(end)]);
set(gca,'Ydir','normal') 
%plot([Ax Bx],[Ay By],'whiteo','MarkerSize',10);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomGLSMZoomk20.eps')


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2SumReg25(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedJReg25.eps')

figure,
hold on, 
contourf(BXmin:DB:BXmax,BYmin:DB:BYmax,J2SumReg25(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2SumReg35(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedJReg35.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2SumReg40(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedJReg40.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2SumReg50(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedJReg50.eps')




figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2SumReg25(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedIReg25.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2SumReg35(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedIReg35.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2SumReg40(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedIReg40.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2SumReg50(:,:,iter_k))
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFixedIReg50.eps')


%% Straight Cracks

load('../TestFF2019/FixedFrq/OutputTemp.mat');

% Les fissures
CentreX = [-2,  2-0.3, 2+0.3,    -2-0.2, -2, -2+0.2,    2-0.4, 2-0.2, 2, 2+0.2, 2+0.4];
CentreY = [ 2,  2    , 2    ,    -2    , -2, -2    ,    -2   ,-2    ,-2,  -2  ,   -2];
r = 0.25;
alpha = pi/2.;
cosinus = r/2*cos(alpha);
sinus = r/2*sin(alpha);
Ax = cosinus + CentreX;  
Ay = sinus + CentreY;
Bx = -cosinus + CentreX;
By = -sinus + CentreY;

Rin = 0.1;
BXmin = -3.125;
BXmax = 3.125;
BYmin = -3.125;
BYmax = 3.125;
DB = Rin/2.; % pas du balayage
[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum)
xlim([-3 3]);
ylim([-3 3]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/StraightFixedJReg2.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2Sum)
xlim([-3 3]);
ylim([-3 3]);
set(gca,'Ydir','normal')
figCracks = plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/StraightFixedIReg2.eps')




%% Single crack demo
Rin = .3;
Folder = 'FixedFrqData/1Crack/';
    load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'I1','J1','I2','J2',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By')
disp([Folder 'FFRin' num2str(round(10*Rin)) '.mat' ' loaded'] );


iter_k = 3;
indic = 3;

CrackColor = Color.gray2;
Crackeps = 3;

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([-0.65 0.55]);
ylim([-0.65 0.55]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedIR3.eps')


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([-0.65 0.55]);
ylim([-0.65 0.55]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedJR3.eps')

%% Single crack demo
Rin = .1;
Folder = 'FixedFrqData/1Crack/';
    load([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],'I1','J1','I2','J2',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By','List_k')
disp([Folder 'FFRin' num2str(round(10*Rin)) '.mat' ' loaded'] );


[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nk = length(List_k)
Nt = numel(BX)
iter_k = 3;
indic = 3;

CrackColor = Color.gray2;
Crackeps = 3;

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedIR1.eps')


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedJR1.eps')

% reg

disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
RegList = sort([1:100 sqrt(2)*(1:100)]);

for RegLevel = 3

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



    for iter_k = 3;

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

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedJR1Reg1.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2Sum(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);

saveanothername([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],J2Sum,['J2SumReg' num2str(RegLevel)]);
saveanothername([Folder 'FFRin' num2str(round(10*Rin)) '.mat'],I2Sum,['I2SumReg' num2str(RegLevel)]);

end

%% Single crack demo
Rin = .01;
Folder = 'FixedFrqData/1Crack/';
    load([Folder 'FFRin01' '.mat'],'I1','J1','I2','J2',...
        'Rin','BXmin','BXmax','BYmin','BYmax','DB','Ax','Ay','Bx','By','List_k')
disp([Folder 'FFRin' num2str(round(10*Rin)) '.mat' ' loaded'] );


[BX,BY]=meshgrid(BXmin:DB:BXmax,BYmin:DB:BYmax);
Nk = length(List_k)
Nt = numel(BX)
iter_k = 3;
indic = 3;

CrackColor = Color.gray2;
Crackeps = 3;

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedIR01.eps')


figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedJR01.eps')

% reg

disp(['~~~~~~~~~~ REGULARIZATION ~~~~~~~~~~']);
RegList = sort([1:100 sqrt(2)*(1:100)]);

for RegLevel = 1:6

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



iter_k = 3
    parfor iter_t = 1:Nt
    
    disp([num2str(iter_t) '/' num2str(Nt)]);

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


toc

%

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,J2Sum(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
%print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/DemoFixedJR1Reg01.eps')

figure,
hold on, 
imagesc(BXmin:DB:BXmax,BYmin:DB:BYmax,I2Sum(:,:,iter_k))
set(gca,'Ydir','normal')
plot([Ax;Bx],[Ay;By],'LineWidth',Crackeps,'Color', CrackColor);
axis square
colormap(jet)
colorbar
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);


saveanothername([Folder 'FFRin01'  '.mat'],J2Sum,['J2SumReg' num2str(RegLevel)]);
saveanothername([Folder 'FFRin01'  '.mat'],I2Sum,['I2SumReg' num2str(RegLevel)]);

end
%% FM Crack

load('FixedFrqData/random2/FFk50');
lambda = 2*pi./k
[nbcapteur,~,Nk] = size(FF);
BXmin = -2.5; BXmax = 2.5;
BYmin = -2.5; BYmax = 2.5;
clc
eta=.035;
Feta = zeros(nbcapteur,nbcapteur,Nk);
Fart = zeros(nbcapteur,nbcapteur,Nk);
HEST = zeros(1,Nk);
CrackColor = Color.gray2;
Crackeps = 3;

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
% Factorisations method
global Feta HEST List_k

iter_k = 1
    ListX = BXmin:.01:BXmax;
    ListY = BYmin:.01:BYmax;
 Mreg = 500

I_FM = CrackFM(ListX,ListY,iter_k,Mreg);
figure,
hold on,
imagesc(ListX,ListY,I_FM)
xlim([-2.5 2.5]);
ylim([-2 2.5]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFMk50.eps')

%zoom

load('FixedFrqData/random2/FFk50');

[nbcapteur,~,Nk] = size(FF);


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


BXminZoom = 0.4;
BXmaxZoom = 1.4;
BYminZoom = -1.5;
BYmaxZoom = -0.5;


ListXZoom = BXminZoom:0.001:BXmaxZoom;
ListYZoom = BYminZoom:0.001:BYmaxZoom;

% Fissures distribution zoom
BXminZoom = ListXZoom(1);
BXmaxZoom =  ListXZoom(end);
BYminZoom = ListYZoom(1);
BYmaxZoom =  ListYZoom(end);
figure,hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
figProbe = plot([BXminZoom;BXmaxZoom;BXmaxZoom;BXminZoom;BXminZoom],[BYminZoom;BYminZoom;BYmaxZoom;BYmaxZoom;BYminZoom],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',25)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/CrackDistributionZoom.eps')




Mreg = 500
I_FM = CrackFM(ListXZoom,ListYZoom,iter_k,Mreg);
figure,
hold on,
imagesc(ListXZoom,ListYZoom,I_FM)
xlim([BXminZoom BXmaxZoom]);
ylim([BYminZoom BYmaxZoom]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFMkZoom50.eps')



%%

load('FixedFrqData/random2/FFk21');
lambda = 2*pi./k
[nbcapteur,~,Nk] = size(FF);
BXmin = -2.5; BXmax = 2.5;
BYmin = -2.5; BYmax = 2.5;

eta=.017;
Feta = zeros(nbcapteur,nbcapteur,Nk);
Fart = zeros(nbcapteur,nbcapteur,Nk);
HEST = zeros(1,Nk);
CrackColor = Color.gray2;
Crackeps = 3;

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
% Factorisations method
global Feta HEST List_k

iter_k = 1
    ListX = BXmin:.01:BXmax;
    ListY = BYmin:.01:BYmax;
 Mreg = 500

I_FM = CrackFM(ListX,ListY,iter_k,Mreg);
figure,
hold on,
imagesc(ListX,ListY,I_FM)
xlim([-2.5 2.5]);
ylim([-2 2.5]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version7/Images/RandomFMk21.eps')


%% FM Higher frq

%% FM Crack

load('FixedFrqData/random2/FFk120');

[nbcapteur,~,Nk] = size(FF);
BXmin = -2.5; BXmax = 2.5;
BYmin = -2.5; BYmax = 2.5;
clc
eta=.035;
Feta = zeros(nbcapteur,nbcapteur,Nk);
Fart = zeros(nbcapteur,nbcapteur,Nk);
HEST = zeros(1,Nk);
CrackColor = Color.gray2;
Crackeps = 3;

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
% Factorisations method
global Feta HEST List_k

iter_k = 1
    ListX = BXmin:.01:BXmax;
    ListY = BYmin:.01:BYmax;
 Mreg = 1000

I_FM = CrackFM(ListX,ListY,iter_k,Mreg);
figure,
hold on,
imagesc(ListX,ListY,I_FM)
xlim([-2.5 2.5]);
ylim([-2 2.5]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFMk120.eps')

%%
%zoom

load('FixedFrqData/random2/FFk120');
lambda = 2*pi./k
[nbcapteur,~,Nk] = size(FF);


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


BXminZoom = 0.4;
BXmaxZoom = 1.4;
BYminZoom = -1.5;
BYmaxZoom = -0.5;


ListXZoom = BXminZoom:0.001:BXmaxZoom;
ListYZoom = BYminZoom:0.001:BYmaxZoom;

% Fissures distribution zoom
BXminZoom = ListXZoom(1);
BXmaxZoom =  ListXZoom(end);
BYminZoom = ListYZoom(1);
BYmaxZoom =  ListYZoom(end);
figure,hold on 
figCracks = plot([Ax;Bx],[Ay;By],'black','LineWidth',Crackeps);
figProbe = plot([BXminZoom;BXmaxZoom;BXmaxZoom;BXminZoom;BXminZoom],[BYminZoom;BYminZoom;BYmaxZoom;BYmaxZoom;BYminZoom],'blue','LineWidth',5);
legend([figCracks(1) figProbe(1)],'Cracks','Probed region','FontSize',25)
xlim([BXmin BXmax]);
ylim([BYmin BYmax]);
set(gca,'Ydir','normal')
axis square
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/CrackDistributionZoom.eps')




Mreg = 500
I_FM = CrackFM(ListXZoom,ListYZoom,iter_k,Mreg);
figure,
hold on,
imagesc(ListXZoom,ListYZoom,I_FM)
xlim([BXminZoom BXmaxZoom]);
ylim([BYminZoom BYmaxZoom]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFMkZoom120.eps')


%% FM Higher frq ++

%% FM Crack

load('FixedFrqData/random2/FFk628');
lambda = 2*pi./k
[nbcapteur,~,Nk] = size(FF);
BXmin = -2.5; BXmax = 2.5;
BYmin = -2.5; BYmax = 2.5;
clc
eta=.035;
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
% Factorisations method
global Feta HEST List_k

iter_k = 1
    ListX = BXmin:.01:BXmax;
    ListY = BYmin:.01:BYmax;
 Mreg = 500

I_FM = CrackFM(ListX,ListY,iter_k,Mreg);
figure,
hold on,
imagesc(ListX,ListY,I_FM)
xlim([-2.5 2.5]);
ylim([-2 2.5]);
set(gca,'Ydir','normal') 
axis square
colormap(jet)
colorbar
print(gcf,'-depsc2','/Users/kev/Downloads/version6/Images/RandomFMk628.eps')