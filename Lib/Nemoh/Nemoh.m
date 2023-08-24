%--------------------------------------------------------------------------------------
%
%    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Contributors list:
%   - A. Babarit
%   - R. Kurnia
%--------------------------------------------------------------------------------------
%
% --> function [Idw,w,A,B,Fe]=Nemoh(projdir,ID_HydrosCal,ID_QTF)
%
% Purpose: Matlab wrapper for calculation of hydrodynamic coefficients using Nemoh
%
% Inputs :
% - projdir     : Project directory path
% - ID_HydrosCal: A switch to compute hydrostatics, 1 computed, 0 Not
% - ID_QTF      : A switch to compute QTF, 1 computed, 0 Not
% Outputs :
% - Idw: A frequency type: 1,2,3= rad/s,Hz,s
% - w  : Vector length(w) of wave frequencies (rad/s)
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of exciation forces (complex
% values)
%
function [Idw,w,A,B,Fe]=Nemoh1(projdir,ID_HydrosCal,ID_QTF)
system(['mkdir ',projdir]);
system(['mkdir ',projdir,filesep,'mesh']);
system(['mkdir ',projdir,filesep,'results']);
if ID_HydrosCal==1
    system(['mkdir ',projdir,filesep,'Mechanics']);
end
if ID_QTF==1
    system(['mkdir ',projdir,filesep,'Motion']);
    system(['mkdir ',projdir,filesep,'results',filesep,'sources']);
end

% Calcul des coefficients hydrodynamiques
fprintf('\n------ Starting NEMOH ----------- \n');
system(['preProc ',projdir]);
if ID_HydrosCal==1
    fprintf('------ computes Hydrostatic ------------- \n');
    system(['hydrosCal1 ',projdir]);
end
fprintf('------ Solving BVPs ------------- \n');
system(['solver ',projdir]);
fprintf('------ Postprocessing results --- \n');
system(['postProc ',projdir]);
%% Lecture des resultats CA CM Fe
%clear Periode A B Famp Fphi Fe;

fid=fopen([projdir,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:nBodies
    for ii=1:4
        ligne=fgetl(fid);
    end
    Ndof=fscanf(fid,'%g',1);
    for j=1:Ndof
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
    Nforc=fscanf(fid,'%g',1);
    for j=1:Nforc
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
ligne=fgetl(fid);
ligne=fscanf(fid,'%g',2);
Idw=ligne(1);
nw=ligne(2);
fclose(fid);
fid=fopen([projdir,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
for c=1:Nforc*nBodies
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+2*Nforc*nBodies);
    w(k)=ligne(1);
    for j=1:Nforc*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end
end
status=fclose(fid);
fid=fopen([projdir,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:Ndof*nBodies
    ligne=fgetl(fid);
end
for i=1:nBodies*Ndof
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
        for j=1:Ndof*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end
        ligne=fgetl(fid);
    end
end
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
Fe=Famp.*(cos(Fphi)+1i*sin(Fphi));
end
