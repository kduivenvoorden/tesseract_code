function [CNOTS_X,CNOTS_Z]  = createCNOTS(N,rs,f)
%determines the CNOTS needed to measure all stabilizers in a tesseract code

%INPUT
% N     size of lattice L = 2^N+2
% rs    boundaries (r = rough) (s = smooth)
% f     possible extension in time dimension

%OUTPUT
% CNOTS_X   a cell of 8 matrices of size nq*2. Each cell represents a set
%           of CNOTS which can be applied in parrallel. They are applied
%           between a qubit and a stabilizer ancilla, which are stored as
%           a tuple, (qubit,stab). Each round a maximum of nq CNOTS are applied.
% CNOTS_Z   Simmilarly as CNOTS_Z

% current lattice (without time direction)
global lattice_s
% length of lattice
L           =2^N+2;
% dimension of lattice
D           = length(rs);
% definition of qubits (in terms of the four corners of a face)
qubit_def   = lattice_s.dCells_all{3,N+1};
%number of qubits in lattice
nq          = size(qubit_def,1);
% paritiy check matrix of X-stabilizers
SX_def = lattice_s.dCells_all{2,N+1};
% paritiy check matrix of Z-stabilizers
if D > 2
    SZ_def  = lattice_s.dCells_all{4,N+1};
else
    SZ_def  = 0;
end


% a string used for storage purposes
switch length(rs)
    case 2, rs_str = sprintf('%d%d',rs(1),rs(2));
    case 3, rs_str = sprintf('%d%d%d',rs(1),rs(2),rs(3));
    case 4, rs_str = sprintf('%d%d%d%d',rs(1),rs(2),rs(3),rs(4));
    case 5, rs_str = sprintf('%d%d%d%d%d',rs(1),rs(2),rs(3),rs(4),rs(5));
end
fname = strcat(sprintf('CNOT/CNOT_%d_%d_',N,f),rs_str,'.mat');


%check if CNOTS are already created...
if exist(fname,'file')==0    
    %... if not, create them
    
    %initialeze the cell of matrices
    CNOTS_X = cell(2*D,1);
    CNOTS_Z = cell(2*D,1);
    
    %the cells contain those CNOTS which can be peformed in parrallel since
    %they are all oriented in a certain direcion
    for i=1:2*D
        % the orientation of cnots depends on i according to the table
        % a = [-1 , -2, -3, -4 , ( ... ) , 4 ,3, 2, 1]
        % the number determining the dimension, the sign the direction
        
        
        %initialeze the matrices
        CNOTS_X{i} = zeros(nq,2);
        CNOTS_Z{i} = zeros(nq,2);
        %1st column: qubit
        %2nd column: stabilizer
        
        %keep track of number of CNOTS, this could be less then nq.
        nx = 0;
        nz = 0;
        
        %starting from the center of a face, and moving in the direction a
        %(uniquely determined by i) one obtains an edge or cube, depending
        %on the orientation of the face. Alternatively, we start with a
        % face f = [v1,v2,v3,v4] and construct the cube 
        % c = [v1+a,v2+a,v3+a,v4+a] = (f+a) in case of Z-stab and check if it is a 
        % ligitimate cube. Or start with an 
        % edge e = [v1,v2] and construct the face
        % f = [v1,v2,v1-a,v2-a] = (e-a) in case of X-stab
        
        % vertices are indices of the form 1+\sum_i (c_i-1)L^(i-1) so a shift in the
        % i direcion is equivalent to adding L^(i-1)
        
        
        if i<=D
            % direction is negative, so consider 
            % the face  f=e+L^.. or 
            % the cube  f=e-L^..  
            X = [SX_def,   SX_def+L^(i-1)];
            Z = [qubit_def, qubit_def-L^(i-1)];
        else
            % direction is negative, so consider 
            % the face  f=e-L^.. or 
            % the cube  f=e+L^.. 
            X = [SX_def,   SX_def-L^(2*D-i)];
            Z = [qubit_def, qubit_def+L^(2*D-i)];
        end
        
        % check for all f=e+a if it is a ligitimate face
        % for all faces
        for j=1:nq
            % which rows of X are equal to such a face
            x = find(all(ismember(X,qubit_def(j,:)),2),1);
            % corresponding row (x) is the stabilizer-def used to connect
            % to this phase
            if ~isempty(x)
                % increase CNOT counter
                nx = nx+1;
                % add pair qubit-stab
                CNOTS_X{i}(nx,:) = [j,x];
            end
        end

        % check for all c=e+a if it is a ligitimate cube
        % for all cubes
        for j=1:size(SZ_def,1)
            % which rows of Z are equal to such a cube
            x = find(all(ismember(Z,SZ_def(j,:)),2),1);
            % corresponding row (x) is the qubit-def used to connect
            % to this phase
            if ~isempty(x)
                % increase CNOT counter
                nz = nz+1;
                % add pair qubit-stab
                CNOTS_Z{i}(nz,:) = [x,j];
            end
        end

        %only keep that part of hte matrix which is actually storing
        %something
        CNOTS_X{i} = CNOTS_X{i}(1:nx,:);
        CNOTS_Z{i} = CNOTS_Z{i}(1:nz,:);
    end
    % save the CNOTS for future use.
    save(fname,'CNOTS_X','CNOTS_Z');
else
    % if CNOTS where already created, just load the corresponding file
    load(fname)
end
end