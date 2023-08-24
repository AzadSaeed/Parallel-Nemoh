% This Script demonstrates parallel computation with Nemoh when 
% one attempts to solve a large number of cases. This code does not change
% the inner-workings of nemoh, rather, it utilizes the capabilities of your
% machine to perform multiple runs in parallel. 
%
% The Matlab wrappers created by Nemoh developers are used in this 
% implementation,, but some of them are slightly modified to account for
% input variations.
%
% The package includes all the required Nemoh files and there is no need to
% download Nemoh, separately.
% 
% The example below is a simple cylinder with variations in radius, only.
% You can easily change the specific nemoh function (i.e. NemohOptions_S )
% to increase the number of cases that run in parallel, or define
% additional cases. 
%
%
% Primary contributor: Dr. Saeed Azad
% Email: saeed.azad@colostate.edu
% https://github.com/AzadSaeed/Parallel-Nemoh

clear; clc; close all;



% Directories:

% Working Directory
pdata.working_dir = pwd;

% Craete a storage directory fo saving data
pdata.Storage_dir = 'Data';
if ~isfolder(pdata.Storage_dir)
    mkdir(pdata.Storage_dir)
end

% Create a plot directory for saving figures
pdata.out_dir = 'Output';
if ~isfolder(pdata.out_dir)
    mkdir(pdata.out_dir)
end


% Nemoh Options - General 
pdata = NemohOptions_G(pdata);


% Nemoh Options - Specific 
pdata = NemohOptions_S(pdata);


% Implement parallel computation
p = gcp();
for idx = 1:pdata.ncases 
  f(idx) = parfeval(p,@(p)ManagePar(p,pdata),1,idx); 
end

% Collect the results as they become available.
NemohResults = cell(1,pdata.ncases);

for idx = 1:pdata.ncases

    % fetchNext blocks until next results are available.
    [completedIdx,value] = fetchNext(f);
    NemohResults{completedIdx} = value;
    fprintf('Got result with index: %d.\n', completedIdx);

end

% Save data in the storage directory
save(strcat(pdata.Storage_dir, filesep,'results'),'NemohResults')


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function pdata = NemohOptions_G(pdata)
% General Nemoh options: these are parameters that are fixed in all of the
% cases. 

% Number of wave energy converters
pdata.nwec = 1;

% Mesh method: a
pdata.MeshMethod = 'General';

% Number of frequency components
pdata.nbfreq     = 10 ;

% Minimum frequency
pdata.wMin     = 0.1;

% Maximum frequency
pdata.wMax     = 2;

% Frequency vector
pdata.w        = linspace(pdata.wMin,pdata.wMax,pdata.nbfreq);

% Water depth m
pdata.depth      = 50; 

% Water density Salt water density kg/m^3 
pdata.rho        = 1025;              

% Gravitational constant
pdata.g          = 9.81;


end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function pdata = NemohOptions_S(pdata)
% Specific Nemoh options: these are parameters that change for every case.

% wave energy converter radius vector
pdata.WEC_radius = linspace(0.5,15,20);

% wave energy converter draft vector
pdata.WEC_draft  = -0.8*ones(1,length(pdata.WEC_radius));

% Total number of cases to run
pdata.ncases     = length(pdata.WEC_radius);

end



