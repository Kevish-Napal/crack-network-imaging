%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtHmxHelmholtz2dDt.m                         |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Martin Averseng             |
%|  ( # )  |   CREATION   : 17.02.2020                                   |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Create FF for crack networks using double     |
%|  `---'  |                layer transpose potential                     |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('/Users/kev/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')

% Parameters
N   = 1e3
tol = 1e-3
typ = 'P1'
gss = 3


%FileOutput = 'FixedFrqData/random2/FFk21.mat';
%List_k = [21]

 
 List_k = [15:0.1:40]

nbcapteur = 100;
Nk = length(List_k)
lambda = 2*pi./List_k

for iterNbCracks = 1:10
    iterNbCracks
FileOutput = ['MultiFrqData/Monotonicity/FF' num2str(iterNbCracks) '.mat'];

% %choose geometry
% CentreX_Full = [-0.5,2.,-1.5,-0.5,  -1., -0.8, -0.75,-1.1,1.5,1.5,1.3,1.6,1.3,1.2];
% CentreY_Full = [ 2, .5, .5, 0.25,  -0.7, -1.1,- 1.7,-2,-1.5,-1.4,-1.6,-1.6,-1.4,-1.5];
% rCrack_Full = [1.5,3,  1., 1, 1, 0.5,0.5,0.75,0.75,0.25,0.25,0.25,0.25,0.5]; 
% alphaCrack_Full = [1.2,0.75,  2/3., 0.15,  0.55, 0.55, 0.3,0.4,0.95,0.9,1,1.05,0.9,0.95];
% alphaCrack_Full = pi*alphaCrack_Full;
% SelectCrack = 1:14;%length(CentreX_Full);
% %select
% CentreX_Full(SelectCrack);
% CentreY_Full(SelectCrack);
% rCrack_Full(SelectCrack);
% alphaCrack_Full(SelectCrack);


% %random cracks
% NbCracks = 40;
% CentreX = [-1+randn(1,NbCracks/4),1+0.5*randn(1,NbCracks/4),-1+0.4*randn(1,NbCracks/4),1+0.2*randn(1,NbCracks/4)];
% CentreY = [1+randn(1,NbCracks/4),1+0.5*randn(1,NbCracks/4),-1+0.4*randn(1,NbCracks/4),-1+0.2*randn(1,NbCracks/4)];
% rCrack = abs(0.2+0.2*randn(1,NbCracks))
% alphaCrack = pi*rand(1,NbCracks);

%vertical cracks
% CentreX = [0];
% CentreY = [0.,0.];
% rCrack = [1.,1.];
% alphaCrack = pi/2.*[1,1];

% or 
% Ax = [0];
% Ay = [0.25];
% Bx = [0];
% By = [-0.25];


% cosinus = rCrack/2.*cos(alphaCrack);
% sinus = rCrack/2.*sin(alphaCrack);
% Ax = cosinus + CentreX;  
% Ay = sinus + CentreY;
% Bx = -cosinus + CentreX;
% By = -sinus + CentreY;

% load existant crack
%load('FixedFrqData/random2/FF.mat','Ax','Ay','Bx','By');
load(['MultiFrqData/Monotonicity/geometry' num2str(iterNbCracks) '.mat']);


% Boundary mesh
% mesh = mshCrack(CentreX,CentreY,rCrack,alphaCrack,N);
mesh = mshCrackBis(Ax,Ay,Bx,By,N);
% Frequency adjusted to maximum esge size
stp  = mesh.stp;
kmax = 1/stp(2)
lambdamin = 2*pi/kmax


% plot(mesh);
% axis square

%

FF = zeros(nbcapteur,nbcapteur,Nk);

% Plane waves direction (to compute FF only)
theta = 2*pi/nbcapteur .* (0:nbcapteur-1)';
nu    = [cos(theta),sin(theta),zeros(size(theta))];




for iter_k = 1:Nk
    tic
    disp([num2str(iter_k), '/', num2str(Nk)]) ;
    k = List_k(iter_k);
    
if (k > kmax)
    warning('Wave number is too high for mesh resolution')
end
f = (k*340)/(2*pi);


%%% PREPARE OPERATOR
%disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Green kernels
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);

% Domain
sigma = dom(mesh,gss);

% Finite elements
u = fem(mesh,typ);
v = fem(mesh,typ);

% Mass matrix
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dnyG(x,y) psi(y) dx dy 
D = (1i/4) .* integral(sigma,sigma,u,dyGxy,ntimes(v),tol);


% Regularization
Dr = -1/(2*pi) .* regularize(sigma,sigma,u,'grady[log(r)]',ntimes(v));

% Operator [-Id/2 + Dt]
LHS = - 0.5*Id + (D + Dr).';

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Sinf = exp(1.i*pi/4.)/sqrt(8.*pi*k) * integral(nu,sigma,Ginf,v);


for itercapteur = 1:nbcapteur


% Incident wave
X0  = [cos(theta(itercapteur)) sin(theta(itercapteur)) 0];
PW = @(X) exp(1i*k*X*X0');
gradxPW{1} = @(X) 1i*k*X0(1) .* PW(X);
gradxPW{2} = @(X) 1i*k*X0(2) .* PW(X);
gradxPW{3} = @(X) 1i*k*X0(3) .* PW(X);




% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
RHS = - integral(sigma,ntimes(u),gradxPW);


%%% SOLVE LINEAR PROBLEM
%disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% LU factorization
[Lh,Uh] = lu(LHS);

% subplot(2,2,3)
% spy(Lh)
% subplot(2,2,4)
% spy(Uh)

% Solve linear system [-Id/2 + Dt] * lambda = dnP0

lambda = Uh \ (Lh \ RHS); % LHS \ RHS;



%%% INFINITE SOLUTION
%disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Finite element radiation  
sol = Sinf * lambda;

FF(:,itercapteur,iter_k) = sol;


end %fin de boucle remplissage du FF
    toc
    
end
save(FileOutput,'FF','k','Ax','Ay','Bx','By','List_k');


disp('~~~~~~~~~~~~~ DONE ~~~~~~~~~~~~~')

end 
%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Radiative mesh (to compute solution)
LSquare = 7;
radiat = mshSquare(2*LSquare*N,[LSquare LSquare]);

% Mesh representation
figure
plot(mesh)
hold on
plot(radiat)
plotNrm(mesh)
axis equal
axis(3.5*[-1 1 -1 1 -1 1])

% Incident wave representation
plot(radiat,real(PW(radiat.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
alpha(0.80)
view(0,90)


% Green kernels
Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1i/4 .* integral(radiat.vtx,sigma,Gxy,v,tol);
toc

% Regularization
tic
Sreg = -1/(2*pi) .* regularize(radiat.vtx,sigma,'[log(r)]',v);
Sdom = Sdom + Sreg;
toc

% Domain solution
Psca = Sdom * lambda;
Pinc = PW(radiat.vtx);
Pdom = Psca + Pinc;

% Incident wave representation
FigUi = figure,hold on 
plot(radiat,real(Pinc))
plot([Ax;Bx],[Ay;By],'red','LineWidth',2);
axis equal
hold off
print(FigUi,'-depsc2','/Users/kev/Downloads/Poster/Ui.eps')
colorbar
title('Incident wave')



% Graphical representation utot
FigUtot = figure,hold on 
plot(radiat,real(Pdom))
plot([Ax;Bx],[Ay;By],'red','LineWidth',2);
axis equal
hold off
print(FigUtot,'-depsc2','/Users/kev/Downloads/Poster/Utot.eps')
title('Total field solution')
colorbar


% Graphical representation us
FigUs = figure,hold on 
plot(radiat,real(Psca))
plot([Ax;Bx],[Ay;By],'red','LineWidth',2);
axis equal
hold off
print(FigUs,'-depsc2','/Users/kev/Downloads/Poster/Us.eps')
title('scattered field')
colorbar


%% Graphical representation us
FigGeometry = figure 
plot([Ax;Bx],[Ay;By],'LineWidth',7,'Color', [.5 .5 .5]);
axis equal
axis([-2.5 3.5 -3 3])
print(FigGeometry,'-depsc2','/Users/kev/Downloads/Poster/Geometry.eps')



disp('~~> Michto gypsilab !')


