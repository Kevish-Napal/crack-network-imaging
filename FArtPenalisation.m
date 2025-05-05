function [gSharp,gTK,aTK,Penalisation] = FArtPenalisation(Rin,t,BC,k,Feta,h,DiscreteDiskCentered)

%global DiscreteDiskCentered Feta HEST List_k


[nbcapteur,~] = size(Feta);
capteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;


% Translation du farfield
ThetaCapteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;
[d,xhat] = meshgrid(ThetaCapteur);
[Xd,Yd] = pol2cart(d,1);
[Xxhat,Yxhat] = pol2cart(xhat,1);
sx = Xd - Xxhat; sy = Yd - Yxhat;
T = @(k,t) exp(1.i*k*(t(1)*sx + t(2)*sy));

FDiskCentered = FGen(k,Rin,[0 0],nbcapteur,BC);
Fbt = T(k,t).*FDiskCentered;
Fart = Feta - Fbt;

DiscreteDisk = DiscreteDiskCentered + t';

% pour tikhonov sur Fart g = Phi_z^\infty:

[Ueta,Seta,Veta]=svd(Fart);             % Ueta et Veta sont des matrices de rotations, Seta une matrice diagonale contenant les valeurs singulières de Feta.
SIGMAN=diag(Seta);                      % SIGMAN est un vecteur contenant la diagonale de Seta.
SQ=SIGMAN.^2;                           % Vecteur contenant les valeurs propres de (Feta)'*Feta.
Us=Ueta';
NbPoints = length(DiscreteDisk(1,:));
options=optimset('Display','off');                                          


aTK = zeros(1,NbPoints);
gSharp = zeros(nbcapteur,NbPoints);
gTK = zeros(nbcapteur,NbPoints);

dx=cos(capteur);dy=sin(capteur);               % coordonnées cartésiennes des capteurs.

%%H = HArtSoftDisk(k,Rin,t,nbcapteur);
%%nH = norm(H);

% Calcul de la penalisation (somme des deux Fsharp)
[Vr,Dr]=eig((Feta+Feta')/2);
[Vi,Di]=eig((Feta-Feta')/(2*sqrt(-1)));
Fsharp=Vr*abs(Dr)*Vr^-1+Vi*abs(Di)*Vi^-1;
NFSharp = norm(Fsharp);

[Vr,Dr]=eig((Fbt+Fbt')/2);
[Vi,Di]=eig((Fbt-Fbt')/(2*sqrt(-1)));
Fbtsharp=Vr*abs(Dr)*Vr^-1+Vi*abs(Di)*Vi^-1;
NFbtSharp = norm(Fbtsharp);

Penalisation = Fsharp + Fbtsharp;
        
    for i = 1:NbPoints 
    %% première itération donnée par la LSM
     
        z = DiscreteDisk(:,i);
              
        rhs=transpose(1/sqrt(8.*pi*k)*exp(1i*pi/4.)*(exp(-sqrt(-1).*k.*(dx.*z(1)+dy.*z(2))))); %/norm(exp(-sqrt(-1).*k.*(dx*z(1)+dy*z(2)))));
        sc=Us*rhs;
        GUP=conj(sc).*sc;
        % Definit la fonction dfun ici pour Morozov
        f = @(x) sum( GUP.*(x^2 - (h^2 * SQ))./((SQ+x).^2));
        
        
        
        gammamin=h*min(SIGMAN);    % borne inf pour la recherche de alpha* de morozov.
        gammamax=h*max(SIGMAN);    % le alpha(delta) de morozov se trouve dans l'intervalle [gammamin,gammamax].
      
        
        [a,fval,exitflag] = fzero(f,[gammamin,gammamax],options); %Faire par dichotomie pour voir (nombre de dichotomie?)
        % a =morozov,     fval = dfun(a),     exitflag < 0 => aucuns zero 

        if a<0||exitflag<0
            disp(['Problem: a= ',num2str(a),' Resetting a=0'])
            disp(['Problem: exitflag= ' num2str(exitflag)])
            a=0;
        end

        aTK(i)=a;
        Valxi=(SIGMAN.*sc)./(SQ+a);
        g = Veta*Valxi;
        gTK(:,i) = g;
        
        
    %% Optimisation avec ||Hg||^2
        
%         alpha = a;
%         FH = @(X) alpha * norm(H*(X(:,1) + 1i*X(:,2)))^2 +...
%             norm(Feta*(X(:,1) + 1i*X(:,2)) - rhs)^2;
%         
%         gStarTemp = fminunc(FH,[real(g),imag(g)]);
%         gStar(:,i) = gStarTemp(:,1) + 1i*gStarTemp(:,2);



        aSharp = a/(NFSharp+NFbtSharp);
        gSharp(:,i) = (aSharp*Penalisation+aSharp*h*eye(nbcapteur)+Fart'*Fart)^-1*(Fart'*rhs);


        
    end
    
    

end
