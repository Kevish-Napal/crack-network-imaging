
global DiscreteDisk % discretization of the centered artificial disk

function [I1,I2,J1,J2] = FixedFrqMethod()
% I, J are two different indicators - their definition are provided in the 
% paper. The number 1 or 2 indicates the method % used to compute I or J 
% (1=Tikhonov, 2=GLSM_FSharp).
I1 = zeros(size(BX));
I2 = zeros(size(BX));
J1 = zeros(size(BX));
J2 = zeros(size(BX));


% Sampling points on the centered artificial background
Nbz = 25; % Nb Approx de points de discretisation interieurs du disque artificiel
[X,Y]=meshgrid(linspace(-Rin,Rin,sqrt(4./pi*Nbz)));
Bool = X.^2+Y.^2 < (Rin).^2;
DiscreteDiskCentered = transpose([X(Bool),Y(Bool)]);
NbPoints = length(DiscreteDiskCentered); % Nb de points de discretisation interieurs du disque artificiel
disp(['Discretization of artificial disk: ' num2str(NbPoints) 'points']);
% figure,plot(X(Bool),Y(Bool),'red+'), 
% hold on, plot(Rin*cos(0:0.1:7),Rin*sin(0:0.1:7),'black'),
% axis square;
% title(['Discretisation interieur du disque artificiel centre de rayon ' num2str(Rin) ', nbPtsInterieur: ' num2str(NbPoints)]);
                                                                                                             



% Translation du farfield 
ThetaCapteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;
[d,xhat] = meshgrid(ThetaCapteur);
[Xd,Yd] = pol2cart(d,1);
[Xxhat,Yxhat] = pol2cart(xhat,1);
sx = Xd - Xxhat; sy = Yd - Yxhat;
T = @(k,t) exp(1.i*k*(t(1)*sx + t(2)*sy));

%Etalonnage
tCalibration = [-5,-5];
FDiskCentered = FGen(k,Rin,[0 0],nbcapteur,BC);
Fart = Feta - T(k,tCalibration).*FDiskCentered;
DiscreteDisk = DiscreteDiskCentered + tCalibration';
[gStar,gTK,aTK] = HPenaltyOptimization(k,HEST,Fart,Rin,tCalibration);
H0gTK = HArtSoftDisk(k,Rin,tCalibration,nbcapteur)*gTK;
H0gStar = HArtSoftDisk(k,Rin,tCalibration,nbcapteur)*gStar;

for iter_t = 1:Nt
    tic
    disp([num2str(iter_t) '/' num2str(Nt)]);
    t = [BX(iter_t),BY(iter_t)];
    
    
    DiscreteDisk = DiscreteDiskCentered + t';
    
      % LSM Fart
        Fart = Feta - T(k,t).*FDiskCentered;% ATTENTION: le calcul ici peut etre alleger il suffit de calculer une seule fois FGen avec t=[0,0]!
        
        
        [gStar,gTK,aTK] = HPenaltyOptimization(k,HEST,Fart,Rin,t);
        
        HgTK = HArtSoftDisk(k,Rin,t,nbcapteur)*gTK;
        HgStar = HArtSoftDisk(k,Rin,t,nbcapteur)*gStar;
        
         
  I1(iter_t) = mean(vecnorm(imag(HgTK)));
  I2(iter_t) = mean(vecnorm(imag(HgStar)));
  
  J1(iter_t) = mean(vecnorm(HgTK - H0gTK));
  J2(iter_t) = mean(vecnorm(HgStar - H0gStar));
 
end

end