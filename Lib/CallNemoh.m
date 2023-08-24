function out = CallNemoh(p,pdata)

% Prescribe specifications to solve
Radius  = pdata.WEC_radius(1,p);
Draft   = pdata.WEC_draft(1,p);
Depth   = pdata.depth;

% Create ID.dat file for Nemoh to run properly
ID(pdata);


% Use MESH
x = 0;
y = 0;
z = 0;
pdata.tx  = [x,y,z];
pdata.CG  = zeros(pdata.nwec,3);


[X,pdata] = CreateMeshcood(Radius,Draft,pdata, pdata.MeshMethod);
Np        = pdata.np*ones(pdata.nwec,1);
Nfobj     = pdata.nfobj*ones(pdata.nwec,1);

[Mass, Inertia, KH, XB, YB, ZB] = NemohMesh(pdata.nwec ,Np,X,pdata.tx,pdata.CG,Nfobj,Depth,pdata.w, 0,0,pdata.out_dir_p);


% Define outputs
pdata.Nemoh.Mass = Mass;
pdata.Nemoh.Inertia = [Inertia(1,4,4),Inertia(1,5,5),Inertia(1,6,6)];
pdata.Nemoh.KH  = KH;
pdata.Nemoh.Buoyancy_cood = [XB,YB,ZB];



% Run Nemoh
NemohtimerVal = tic;

% Using Nemoh 3 - not working yet
%[Idw, w, A,B,Fe]      = Nemoh_3(pdata.Nemoh.projdir,pdata.Nemoh.HydrosCal, pdata.Nemoh.QTF);

% Using the older version of Nemoh
[Idw,w,A,B,Fe]      = Nemoh(pdata.out_dir_p,0,0);
Nemohtime           = toc(NemohtimerVal);

pdata.Results.cmpt_Time = Nemohtime;
pdata.Results.A      = A;
pdata.Results.B      = B;
pdata.Results.Fe     = Fe;
pdata.Results.Idw    = Idw;
pdata.Results.w      = w;
pdata.Results.Radius = Radius;
pdata.Results.Draft  = Draft;
pdata.Results.Depth  = Depth;

% declaire output
out = pdata.Results;

end



%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


function ID(pdata)
% Create ID.dat file for Nemoh to run properly
fid = fopen('ID.dat','w');
[~,name,~] = fileparts(pdata.Storage_dir);
fprintf(fid,['% g \n', name,' \n'], length(name));
fclose(fid);
end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function [X,varargout] = CreateMeshcood(radius,draft, pdata, MeshMethod)

switch upper(MeshMethod)

    case 'AXISYMMETRIC'

    % Cylinder
    n              = 3;                     % Number of points of discretization required to define the mesh (3 when axis symmetric (revolution can create the whole mesh))
    r              = [radius radius 0];     % Radial coordinate in the polar system
    z              = [0 draft draft];       % Array of vertical coordinates
    pdata.ntheta   = 60;                    % Number of angular discretization points
    pdata.np       = 504;                   % Number of panels
    pdata.zg       = pdata.CG(1,3);         % Vertical position of gravity center
    

    X = r;
    if nargout >1
        varargout  = {z,n,pdata};
    end

    case 'GENERAL'

    % Cylinder
    pdata.ntheta   = 60;                               % Number of panels on the side
    theta          = 0:pi/(pdata.ntheta-2):pi;         % Create radial panels

    pdata.nbottom  = 90;                               % Number of panels at the bottom of the cylinder
    pdata.nfobj    = 900;                              % Target number of panels for Aquaplus mesh
    pdata.np       = pdata.ntheta-2+pdata.nbottom-1;

    xo = kron(cos(theta),radius);
    yo = kron(sin(theta),radius);

    % Define the coordinates of the radial panels
    x  = zeros(pdata.nwec ,length(xo));
    y  = zeros(pdata.nwec ,length(yo));

    for i = 1:pdata.nwec
        x(i,:) = xo + pdata.CG (i,1);          % Create x coordinate
        y(i,:) = yo + pdata.CG (i,2);          % Create y coordinates
    end
    
    X  = zeros(pdata.nwec,pdata.np,4,3);
    
    for j = 1:pdata.nwec
        for i = 1:pdata.ntheta-2
            X(j,i,:,:) = [x(j,i), y(j,i),draft; x(j,i+1), y(j,i+1), draft;...
                x(j,i+1), y(j,i+1), 0; x(j,i), y(j,i), 0];
        end

        for i = 1:pdata.nbottom
            xb = linspace(xo(end), xo(1), pdata.nbottom);
            yb = sqrt(radius^2 - xb.^2);
            Xb = xb + pdata.CG (j,1);
            Yb = yb + pdata.CG (j,2);
        end

        zb = draft.*ones(1,pdata.nbottom);

        for i = 1:pdata.nbottom-1
            X(j,pdata.ntheta-2+i,:,:) = [Xb(i), Yb(1), zb(i); Xb(i),...
                Yb(i), zb(i); Xb(i+1), Yb(i+1), zb(i); Xb(i+1), Yb(1),...
                zb(i)];
        end
    end

    r = X;
    if nargout >1
        varargout = {pdata};
    end

end


end