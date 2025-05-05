function [I_LSM,I_GLSM] = SamplingMethods(ListX,ListY,iter_k)

global Feta HEST List_k

[nbcapteur,~,~] = size(Feta);
capteur = 0:2*pi/nbcapteur:2*pi-2*pi/nbcapteur;

h = HEST(iter_k);
k = List_k(iter_k);

[ZX,ZY] = meshgrid(ListX,ListY);

% pour tikhonov sur Feta g = Phi_z^\infty:

[Ueta,Seta,Veta]=svd(Feta(:,:,iter_k));             % Ueta et Veta sont des matrices de rotations, Seta une matrice diagonale contenant les valeurs singulières de Feta.
SIGMAN=diag(Seta);                      % SIGMAN est un vecteur contenant la diagonale de Seta.
SQ=SIGMAN.^2;                           % Vecteur contenant les valeurs propres de (Feta)'*Feta.
Us=Ueta';
NbPoints = numel(ZX);
options=optimset('Display','off');                                          


I_LSM = zeros(size(ZX));
I_GLSM = zeros(size(ZX));

dx=cos(capteur);dy=sin(capteur);               % coordonnées cartésiennes des capteurs.

%%H = HArtSoftDisk(k,Rin,t,nbcapteur);
%%nH = norm(H);

% Calcul de Fsharp
[Vr,Dr]=eig((Feta(:,:,iter_k)+Feta(:,:,iter_k)')/2);
[Vi,Di]=eig((Feta(:,:,iter_k)-Feta(:,:,iter_k)')/(2*sqrt(-1)));
Fsharp=Vr*abs(Dr)*Vr^-1+Vi*abs(Di)*Vi^-1;
NFSharp = norm(Fsharp);

        
    parfor i = 1:NbPoints 
    %% première itération donnée par la LSM
     %disp([num2str(i) '/' num2str(NbPoints)]);
        z = [ZX(i) ZY(i)];
              
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
        g = Veta*Valxi;
        
        I_LSM(i) = abs(dot(Feta(:,:,iter_k)*g,g)) + h*dot(g,g);
    %% Optimisation avec ||Hg||^2
        
%         alpha = a;
%         FH = @(X) alpha * norm(H*(X(:,1) + 1i*X(:,2)))^2 +...
%             norm(Feta*(X(:,1) + 1i*X(:,2)) - rhs)^2;
%         
%         gStarTemp = fminunc(FH,[real(g),imag(g)]);
%         gStar(:,i) = gStarTemp(:,1) + 1i*gStarTemp(:,2);



        aSharp = a/NFSharp;
        gSharp = (aSharp*Fsharp+aSharp*h*eye(nbcapteur)+Feta(:,:,iter_k)'*Feta(:,:,iter_k))^-1*(Feta(:,:,iter_k)'*rhs);

        I_GLSM(i) = abs(dot(Fsharp*gSharp,gSharp)) + h*dot(gSharp,gSharp);
        
    end
    
    

end
