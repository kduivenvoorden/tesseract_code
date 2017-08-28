function lattice = createMapsCells(N,rs,f)
%
%creates all edges, faces, cubes, and hypercubes, and their incidence
%matrix and stores it as a struct lattice
%
% using
% createCells.m
% changeBoundary.m
% boundaryOperator.m
%
% N         system size, L = 2^N+2
% rs        rough/smooth, length D vector containing info on which boundaries are
%           rough and which ones are smooth
% f          a parameter for posssibly exteding time direction
%storestyle either global or local


%information is stored not only for system size L(N) = 2^N+2, for also for
%L(N-1), L(N-2) .... , L(0). This information is used in renormalization
%steps. All variables with the extension '_all' are cell arrays such that
% X_all{1}     info of X for L(0) = 2
% X_all{2}     info of X for L(1) = 3
...
% X_all{N+1}   info of X for L(N) = L^N+2

%dCells_all{i,j} contains info of the (i-1)th cell for system size L(j-1)
%SQZ     Z-stabilizers, a binary array of size (nQ,nS) where nS is the
%        number of Z-stabilizers and nQ is the number of qubits
%SQX     X-stabilizers, a binary array of size (nS,nQ) where nS is the
%        number of X-stabilizers and nQ is the number of qubits
%VSX     boundary map of X-stabilizers, binary array of size (nV,nS) where nS is the
%        number of X-stabilizers and nV is the number of vertices

% default values of f, when not given
if nargin ==2
    f = 1;
end

% Dimension of hypercube
D = length(rs);
%initialization of all cell-arrays
dCells_all = cell(D+1,N+1);
SQX_all = cell(N+1,1);
SQZ_all = cell(N+1,1);
VSX_all = cell(N+1,1);

%all cell_arrays are stored in a .mat file after creation. The creation of
%these cell arrays can take up to hours for (N=3 L=10), so it is crucial
%that this information is pre-calculated only once and not at each
%monte-carlo trial


% a string used for storage purposes
switch length(rs)
    case 2, rs_str = sprintf('%d%d',rs(1),rs(2));
    case 3, rs_str = sprintf('%d%d%d',rs(1),rs(2),rs(3));
    case 4, rs_str = sprintf('%d%d%d%d',rs(1),rs(2),rs(3),rs(4));
    case 5, rs_str = sprintf('%d%d%d%d%d',rs(1),rs(2),rs(3),rs(4),rs(5));
end

%for all system sizes
for iN=0:N
    %file name for storing cell definitions
    fname_C = strcat(sprintf('Cells/dCells_%d_%d_',iN,f),rs_str,'.mat');
    %file name for storing maps (stabilizer definitions)
    fname_M = strcat(sprintf('Maps/maps_%d_%d_',iN,f),rs_str,'.mat');
    
    L=2^iN+2;
    %check if cell definitions is already created...
    if exist(fname_C,'file')==0
        %... otherwise, create  cell definitions
        dCells = cell(D+1,1);
        for i=0:D
            % create cell definitions, 
            % all L-(post) boundaries are rough
            % all 1-(pre ) boundaries are smooth 
            C = createCells(D,i,L,f);
            
            %change boundaries
            for j=1:D
                if rs(j)==0
                    %smoothen
                    C = changeBoundary(C,L,D,j,'post','S',f);
                end
                if rs(j)==1
                    %roughen
                    C = changeBoundary(C,L,D,j,'pre','R',f);
                end
            end
            dCells{i+1} = C;   
        end
        %save for further usage
        save(fname_C,'dCells');
     
    end
    
    %check is stabilizer definitions are already created
    if exist(fname_M,'file')==0
        %if not create them
        load(fname_C);
        
        %create maps
        VSX = boundaryOperator(dCells{1},dCells{2});
        SQX = boundaryOperator(dCells{2},dCells{3});
        
        %only define Z-stabilizer if the dimension is big enough (at least
        %3)
        if D > 2
            SQZ = boundaryOperator(dCells{3},dCells{4});
        else 
            SQZ = zeros(size(dCells{3},1),0);
        end
        save(fname_M,'SQX','SQZ','VSX');
    end
    
    % all cell definitions and maps (for this specific system size) are 
    % either created and saved, or were alread created. They are now reloaded
    % and stored in a local variable. (This is inefficient only when
    % creating - saving - reloading the definitions)
    load(fname_C);
    dCells_all(:,iN+1) = dCells;
    load(fname_M);
    SQX_all{iN+1} = SQX;
    SQZ_all{iN+1} = SQZ;
    VSX_all{iN+1} = VSX;
    
end


% definitions are also stored in a struct
lattice.SQX_all  = SQX_all;
lattice.dCells_all = dCells_all;
lattice.SQZ_all = SQZ_all;
lattice.VSX_all = VSX_all;

end

