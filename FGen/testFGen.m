
clear all 
FRealTemp = dlmread(['FRealD.dat']);
FImageTemp = dlmread(['FImageD.dat']);

k = FRealTemp(1,1); Rin = FRealTemp(1,2); ax = FRealTemp(1,3); ay = FRealTemp(1,4);
a = [ax ay];
BC = [1 0];
FF_Real = FRealTemp(2:end,:);
FF_Image = FImageTemp(2:end,:);
FF = FF_Real + sqrt(-1.0)*FF_Image;
nbcapteur = size(FF,1);

PercentageError.Dirichlet = mean(mean(abs(FF-FGen(k,Rin,a,nbcapteur,BC))./abs(FF)))*100;


% Neumann
FRealTemp = dlmread(['FRealN.dat']);
FImageTemp = dlmread(['FImageN.dat']);

k = FRealTemp(1,1); Rin = FRealTemp(1,2); ax = FRealTemp(1,3); ay = FRealTemp(1,4);
a = [ax ay];
BC = [0 1];
FF_Real = FRealTemp(2:end,:);
FF_Image = FImageTemp(2:end,:);
FF = FF_Real + sqrt(-1.0)*FF_Image;
nbcapteur = size(FF,1);

PercentageError.Neumann = mean(mean(abs(FF-FGen(k,Rin,a,nbcapteur,BC))./abs(FF)))*100;


% Impedance
FRealTemp = dlmread(['FRealI.dat']);
FImageTemp = dlmread(['FImageI.dat']);

k = FRealTemp(1,1); Rin = FRealTemp(1,2); ax = FRealTemp(1,3); ay = FRealTemp(1,4);
a = [ax ay];
BC = [1i 1];
FF_Real = FRealTemp(2:end,:);
FF_Image = FImageTemp(2:end,:);
FF = FF_Real + sqrt(-1.0)*FF_Image;
nbcapteur = size(FF,1);

PercentageError.Impedance = mean(mean(abs(FF-FGen(k,Rin,a,nbcapteur,BC))./abs(FF)))*100;
disp('PercentageError:');
disp(PercentageError);
