function [result] = localDecoder(p,N,f,type)
global rs_global visualizes lattice_st lattice_s

% performs a single run of error correction

%INPUT 
%p          probability of error
%N          system size L=2^N+2
%rs         boundaries (r = rough) (s = smooth)
%f          possible extension of time Lt = f(L-2)+2
%type.noise_model   error model, is either fenomological ('fem') or gated based ('gb')
%type.meas_noise    boolean: false if measurements are perfect (no time direction) true otherwise           

%OUTPUT
%result either 1 (failure) or 0 (succes)

%GLOBAL
%rs_global  essentially equal to rs, however, rs can be changed locally
%           (for removing time direcion)
%visualizes boolean, when true the program outputs some graphics for
%           debuggin purposes
%SQX_all, dCells_all
%           defintion of lattice (cells) and stabilizers


%% Preperation
if nargin == 0
    %remove all figures
    close all
    
    %All other .m files are in this folder
    addpath(strcat(pwd,'/Decoder'))
    
    %if it is run as a script (for debugging) the parameters can be set
    %manually
    visualizes = true;
    N=2;
    p = 0.001;
    rs_global = [0;1;1];
    %last dimension is time (and hence smooth)
    type.noise_model = 'fem';
    type.meas_noise = false;
    f=1;
    
    %create the lattice and stabilizer definitions
    lattice_st = createMapsCells(N,rs_global,f);
    lattice_s = createMapsCells(N,rs_global(1:end-1),f);
    
    %------- data for generating specific pictures
    load prob_data_pictures
    syndrome = mod(lattice_st.SQX_all{N+1}*fq,2);
    type.noise_model = '';
    %-------
end

%boundaries of system
rs          = rs_global;
%dimension of system
D           = length(rs);
%spatial dimension of system
D_spatial   = D-type.meas_noise;
%defintion of stabilizers
SQX         = lattice_st.SQX_all{N+1};
%number of qubits
nq          = size(SQX,2);
%qubit weight
wq          = log(p/(1-p))*ones(nq,1);

%create errors at random
switch type.noise_model
    case 'fem'
        %flip qubits,
        fq = rand(nq,1)<p;        
        syndrome = mod(SQX*fq,2);
    case 'gb'
        % weight estimation
        % this function was used to gain info on the estimating the propability
        % distribution of the gated based errors by an independent X/Z
        % model and use the weight 
        
        %wq = weight_estimation(N,f,p,false);
        [CNOTS_X,CNOTS_Z] = createCNOTS(N,rs(1:end-1),f);
        [syndrome, readout_Q_end,~] = createGateBasedErrors(N,f,p,CNOTS_X,CNOTS_Z);  
end
%store the syndrome in a seperate variable, to check later if one has
%mapped back into the code space (new syndrome + old syndrome = 0)
syndrome_original = syndrome;

%% Correction

% N gives the level of renormalization. Errors are removed, first for iN=N
% than for iN=N-1, iN=N-2. Each step takes as input the remaining syndrome of the previous
% output. The qubits that are flipped as a correction are stored.
correction_All = cell(N+1,1);
for iN=N:-1:0
    [syndrome,wq,correction] = renormalization_step(syndrome,wq,iN,f);
    correction_All{iN+1} = correction;
end

% At each level, "re-scaled" qubits can be flipped, a qubit on the iN-th level refers
% to 4 qubits on the (iN+1)-th. The flips on the coarse grained lattices
% are translate to the origal lattice
sc = zeros(size(correction_All{1}));
for iN=0:N-1
    sc = scale_correction(correction_All{iN+1}+sc,iN,f);
end
full_correction = sc+correction_All{N+1,1};


% The should be no syndrome remaining (project back into the code space
if sum(mod(SQX*full_correction+syndrome_original,2))>.5
    fprintf('decoding not succesfull\n')
end

%Prints the flipped qubits (due to errors) and the corrected qubits (due to
%decoding) to check the perfomance of the decoding
if visualizes 
    visualize1(N,full_correction,fq,D,f)
end

%% Analysis of correction

% we will be projecting out the time direction. This will results in a
% lattice of 1 dimension less. The corresponding stabilizers are used to
% check if the logicals are  created correctly. 
qubit_def_l = lattice_s.dCells_all{3,N+1};
SQX = lattice_s.SQX_all{N+1};
QSZ = lattice_s.SQZ_all{N+1};


% create logicals
% Lz is always planar
% Both Lx and Lz is a vector containing the indices of vertices. For Lz it is all vertices
% lying in a certain plain. For Lx, it is all vertices
% lying in a certain plain (4 spatial dimensions) line (D_spatial=3) or merely a single point

L=2^N+2;
if D_spatial ==2
    Lx = 1;
    Lz = 1:max(max(qubit_def_l));
end
if D_spatial ==3
    Lx = my_sub2ind(L, reshape([(1:L)',ones(L,2)],L,1,3),3,1);
    Lz = my_sub2ind(L, reshape([ones(L^2,1),repmat((1:L)',[L,1]),kron((1:L)',ones(L,1))],L^2,1,3),3,1);
end
if D_spatial ==4
     Lx = my_sub2ind(L, reshape([repmat((1:L)',[L,1]),kron((1:L)',ones(L,1)),ones(L^2,2)],L^2,1,4),4,1);
     Lz = my_sub2ind(L, reshape([ones(L^2,2),repmat((1:L)',[L,1]),kron((1:L)',ones(L,1))],L^2,1,4),4,1);
end

% The actual logical will be created by taking those edges which have
% exactly overlap 1 with this plane/line/point (Lx)
Lx = sum(ismember(qubit_def_l,Lx),2)==1;
% or taking all edges in that plane (Lz)
Lz = all(ismember(qubit_def_l,Lz),2);
% Test wether Lz and Lx commute with all stabilizers (and anticommute)
if nargin == 0 && strcmp(type.noise_model,'fem')
    [sum(mod(Lx'*QSZ,2)), sum(mod(SQX*Lz,2)), mod(sum(Lx & Lz),2),sum(fq)]
end

switch type.noise_model
    case 'fem'
        flipped_qubits = full_correction+fq;
    case 'gb'
        flipped_qubits = full_correction;
end

if type.meas_noise
    %project the time direction

    %definition of all faces
    qubit_def = lattice_st.dCells_all{3,N+1};
    faces = my_ind2sub(L,qubit_def);
    %those faces who is perpendicular to the time direction
    isqubit = faces(:,1,D)==faces(:,4,D);

    %remove the time direction of the definition the flipped qubits
    project_qubits = my_ind2sub(L,qubit_def(mod(flipped_qubits,2)==1 & isqubit,:));
    project_qubits = my_sub2ind(L,project_qubits(:,:,(1:end-1)),D-1,1);

    %for each real qubit
    flipped_qubits =zeros(size(qubit_def_l,1),1);
    for i=1:size(project_qubits,1)
        %check how many qubit definitions are fully in projected flipped qubits
        flipped_qubits = flipped_qubits + (sum(ismember(qubit_def_l,project_qubits(i,:)),2)==4);
    end
end

switch type.noise_model
    case 'fem'
        result = mod(Lx'*flipped_qubits,2);
    case 'gb'
        result = mod(Lx'*(flipped_qubits+readout_Q_end),2);
end


end

%% visualization
function visualize1(N,full_correction,fq,D,f)
global lattice_st
L = 2^N+2;
%qubits are surfaces
qubit_def = lattice_st.dCells_all{3,N+1};
if D ==4
    c = mod(full_correction,2)==1; 
    fq = mod(fq,2);
    cf = fq==1 & c==1;

    [x1,y1,z1,u1] = ind2sub([L,L,L,L],qubit_def(fq-cf==1,:));
    [x2,y2,z2,u2] = ind2sub([L,L,L,L],qubit_def(c-cf==1,:));
    [x3,y3,z3,u3] = ind2sub([L,L,L,L],qubit_def(cf,:));    

    for t = 1:L-1

        intializeVisualization(L,f);

        % blue: qubits that are flipped and not corrected
        pt = u1(:,1) == t & u1(:,4)==t;
        plot_e = fill3(x1(pt,[1,2,4,3])'-1,y1(pt,[1,2,4,3])'-1,z1(pt,[1,2,4,3])'-1,[0,0,1]);
        alpha(plot_e,.2);
        % red: qubits that are (wrongly) corrected
        pt = u2(:,1) == t & u2(:,4)==t;
        plot_c = fill3(x2(pt,[1,2,4,3])'-1,y2(pt,[1,2,4,3])'-1,z2(pt,[1,2,4,3])'-1,[1,0,0]);
        alpha(plot_c,.2);
        % green: qubits that are flipped and corrected
        pt = u3(:,1) == t & u3(:,4)==t;
        plot_e = fill3(x3(pt,[1,2,4,3])'-1,y3(pt,[1,2,4,3])'-1,z3(pt,[1,2,4,3])'-1,[0,1,0]);
        alpha(plot_e,.2);
    end
else
    intializeVisualization(L,f);
    c = mod(full_correction,2)==1; 
    fq = mod(fq,2);
    cf = fq==1 & c==1;
    % blue: qubits that are flipped and not corrected
    [x,y,z] = ind2sub([L,L,L],qubit_def(fq-cf==1,:));
    plot_e = fill3(x(:,[1,2,4,3])'-1,y(:,[1,2,4,3])'-1,z(:,[1,2,4,3])'-1,[0,0,1]);
    alpha(plot_e,.2);
    % red: qubits that are (wrongly) corrected
    [x,y,z] = ind2sub([L,L,L],qubit_def(c-cf==1,:));
    plot_c = fill3(x(:,[1,2,4,3])'-1,y(:,[1,2,4,3])'-1,z(:,[1,2,4,3])'-1,[1,0,0]);
    alpha(plot_c,.2);
    % green: qubits that are flipped and corrected
    [x,y,z] = ind2sub([L,L,L],qubit_def(cf,:));
    plot_e = fill3(x(:,[1,2,4,3])'-1,y(:,[1,2,4,3])'-1,z(:,[1,2,4,3])'-1,[0,1,0]);
    alpha(plot_e,.2);
    saveas(gcf,'conclusions.pdf')
end
end

