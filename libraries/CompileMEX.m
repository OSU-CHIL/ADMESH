temp = dir('libraries/admesh library/**/*.c');

PWD = pwd;

for i = 1 : length(temp)

    folder = temp(i).folder;
    name = temp(i).name;

    foldername = [folder,'/',name];

    foldername2 = strrep(foldername,'.c','.mex*');
    
    exist(foldername2,'file')

    if ~exist(foldername2,'file')
        cd(folder);
        mex(foldername);
    end
end

cd(PWD);