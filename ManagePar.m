function out = ManagePar(p,pdata)
% This function manages the parallel implementation of NEMOH by making
% copies of Nemoh source files in local repositories to avoid conflict
% between workers.
%
%
% Primary contributor: Dr. Saeed Azad
% Email: saeed.azad@colostate.edu
% https://github.com/AzadSaeed/Parallel-Nemoh



% Source files
Nemfile1 = 'CallNemoh.m';
srcfile1 = strcat(pdata.working_dir,filesep,'Lib',filesep,Nemfile1);
Nemfile2 = 'FindingNemoh.m';
srcfile2 = strcat(pdata.working_dir,filesep,'Lib',filesep,Nemfile2);

% Destination folder
pdata.out_dir_p = strcat(pdata.working_dir,filesep,'Output',filesep,'Output_',char(string(p)));
if ~isfolder(pdata.out_dir_p)
    mkdir(pdata.out_dir_p)
end

% Copy files to their corresponding folders only if they do not exist
if ~isfile(strcat(pdata.out_dir_p,filesep,Nemfile1))
    copyfile (srcfile1,pdata.out_dir_p)
end

% Copy files to their corresponding folders only if they do not exist
if ~isfile(strcat(pdata.out_dir_p,filesep,Nemfile2))
    copyfile (srcfile2,pdata.out_dir_p)
end

% Copy Nemoh into a local directory
Libpath = strcat(pdata.working_dir, filesep, 'Lib', filesep, 'Nemoh');
Dest = strcat(pdata.out_dir_p,filesep,'Nem_lib');
if ~exist(Dest, 'dir')
    copyfile(Libpath,Dest)
end

% Add all contents to the path in the local directory
addpath(genpath(pdata.out_dir_p));

% Ensure Nemoh is in the path
assert(FindingNemoh(0, true))

% Run Nemoh
out = CallNemoh(p,pdata);


end