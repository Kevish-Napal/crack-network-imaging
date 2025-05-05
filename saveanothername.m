function [] = saveanothername(StringSaveFile,data,StringNewName)
%  StringSaveFile = path/DataFile.mat
% data is saved in the file DataFile.mat
% under the name StringNewName

    S.(StringNewName) = data;

    if isfile(StringSaveFile)
         save(StringSaveFile, '-struct', 'S','-append'); 
    else
         save(StringSaveFile, '-struct', 'S'); 
    end


end