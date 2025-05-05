function ScalaireTemp = FArtPenalisationMulti(Rin,BC,k,Feta,h,DiscreteDiskCentered,BX,BY)


Nt = numel(BX);
Scalaire1Temp = zeros(1,Nt);
Scalaire2Temp = zeros(1,Nt);
Scalaire3Temp = zeros(1,Nt);

[nbcapteur,~] = size(Feta);
capteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;


% Translation du farfield precomputations
ThetaCapteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;
[d,xhat] = meshgrid(ThetaCapteur);
[Xd,Yd] = pol2cart(d,1);
[Xxhat,Yxhat] = pol2cart(xhat,1);
sx = Xd - Xxhat; sy = Yd - Yxhat;

FDiskCentered = FGen(k,Rin,[0 0],nbcapteur,BC);

NbPoints = length(DiscreteDiskCentered(1,:));
options=optimset('Display','off'); 

gSharp = zeros(nbcapteur,NbPoints);
gTK = zeros(nbcapteur,NbPoints);

dx=cos(capteur);dy=sin(capteur);               % coordonnées cartésiennes des capteurs.

for i = 1:NbPoints 
parfor iter_t = 1:Nt

t = [BX(iter_t),BY(iter_t)];



% Translation du farfield
T = @(k,t) exp(1.i*k*(t(1)*sx + t(2)*sy));


Fbt = T(k,t).*FDiskCentered;
Fart = Feta - Fbt;

DiscreteDisk = DiscreteDiskCentered + t';

% pour tikhonov sur Fart g = Phi_z^\infty:

[Ueta,Seta,Veta]=svd(Fart);             % Ueta et Veta sont des matrices de rotations, Seta une matrice diagonale contenant les valeurs singulières de Feta.
SIGMAN=diag(Seta);                      % SIGMAN est un vecteur contenant la diagonale de Seta.
SQ=SIGMAN.^2;                           % Vecteur contenant les valeurs propres de (Feta)'*Feta.
Us=Ueta';
                                         
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

        
        Valxi=(SIGMAN.*sc)./(SQ+a);
        gTK(:,i) = Veta*Valxi;
        
        
        
    %% Optimisation avec ||Hg||^2
        
%         alpha = a;
%         FH = @(X) alpha * norm(H*(X(:,1) + 1i*X(:,2)))^2 +...
%             norm(Feta*(X(:,1) + 1i*X(:,2)) - rhs)^2;
%         
%         gStarTemp = fminunc(FH,[real(g),imag(g)]);
%         gStar(:,i) = gStarTemp(:,1) + 1i*gStarTemp(:,2);



        aSharp = a/(NFSharp+NFbtSharp);
        gStar = (aSharp*Penalisation+aSharp*h*eye(nbcapteur)+Fart'*Fart)^-1*(Fart'*rhs);


         Scalaire1Temp(iter_t) = norm(norm(HArtSoftDisk(k,Rin,t,nbcapteur)*gTK));
        Scalaire2Temp(iter_t) = norm(norm(HArtSoftDisk(k,Rin,t,nbcapteur)*gStar));
        Scalaire3Temp(iter_t) = norm(dot(Penalisation*gStar,gStar));
    end
    
       
    
end

    ScalaireTemp(:,1) = Scalaire1Temp;
    ScalaireTemp(:,2) = Scalaire2Temp;
    ScalaireTemp(:,3) = Scalaire3Temp;
end
