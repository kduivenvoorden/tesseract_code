function B = boundaryOperator(Cells1,Cells2)

% determines the boundary operator \del_d
% by looking at over laps between 
% (d-1)-cells and  d-cells

% Cells1 and Cells2 are matrices containing the indices of the 
% vertices of different cells, for details, check 'createCells.m'


% the output, B will be a [n, m] (binary) matrix whose i-j-th component
% is 1 iff the i-th d1-cells is contained in the j-th d2-cell.
% n = nCells1, number of d1-cells
% m = nCells2, number of d2-cells
[nCells1,s1] = size(Cells1);
[nCells2,s2] = size(Cells2);
%s = 2^d, number of vertices per d-cell

% d1 does not necessarily have to be d2-1, but should at least be smaller
% than d2, otherwise d1-cells can not be contained in d2-cells.
if s2<=s1
    fprintf('Warning: wrong use of script "boundaryOperator"')
end


% The output is stored in as a sparse matrix. 
% To estimate the number of non-zero entries we use that each 
% d-cell contains 2d (d-1) cells.
% note that if d1 ~= d2-1, this upperbound is wrong and will lead to
% inefficient memory allocation.

upperbound = nCells2*2*log2(s2);

% the row and collumn indices of the sparse matrix B
r=zeros(upperbound,1);
c=zeros(upperbound,1);

%this script could take long if n2 is large (>10^4). In this case it ouputs its
%progress
progress=1;
p=100;


%k0 keeps track of how many non-zero entries of B have been found
k0=0;

%per d2 cell ... 
for j=1:nCells2
    % ...check which d1 cell are contained in it...
    Q = all(ismember(Cells1,Cells2(j,:)),2);
    k = sum(Q);
    % ... and store their row index (which is the row index for B)
    r(k0+1:k0+k) = find(Q);
    % The relavent column index for B is just the corresponding d2 cell,
    % ie. j
    c(k0+1:k0+k) = j*ones(k,1);
    k0= k0+k;
    
    %output progress
    if nCells2>10^4
        if j>nCells2/p*progress
            fprintf('%dp',100/p*progress)
            progress=progress+1;
        end
    end
end

% use only the relavent entries of r and c. (The upperbound is generally
% not tight)
r = r(1:k0);
c = c(1:k0);

%create the sparse matrix
B = sparse(r,c,1,nCells1,nCells2,size(c,1));