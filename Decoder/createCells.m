function dCells = createCells(D,d,L,f)
%
% creates d-Cells in a D-dim hypercube of size L, and (possibly extended) time

% d-cells are determined by their vertices
% the output dCells is a [n,s] matrix whose i-j component is the index 
% of the j-th vertix belonging to the i-th cell.
% n = number of cells
% s = 2^d

% for example, if Cells1 are 1-cells (edges) we have that 
% Cells1 = 
%   [v11, v12;
%    v21, v22;
%    ...
%    vn1, vn2]
% where v_i1 and v_i2 are the (indices of the) 2 vertices 
% which are joined by edge i

% D     dimension of the hypercube
% d     dimension of cells
% L     size of hypercube
% f     possibly extending time, the time direction will have size 
%       L(f) = f*(L-2)+2

% if f is not given, set it to 1.
if nargin == 3
    f = 1;
end

% linear indices of vertices run from 1  --- L^D
% for example, D=2, L=3
% 1 4 7
% 2 5 8
% 3 6 9

% d- cells are created by 
% - taking a starting vertex (unit cell) and
% - taking d out of D directions
%   denote them with i_k with k \in [1..D] 

% Let (x1,x2,... xD) be the subscript of the starting vertix. Let y be the 
% corresponding index. The other
% 2^d-1 vertices have linear subscripts and indices
% (x_1,...(x_i1)+1, ... x_i2  , .... x_in  ,.... x_D) -> y + L^(i1-1)
% (x_1,...x_i1  , ...   x_i2+1, .... x_in  ,.... x_D) -> y + L^(i2-1)
% ...
% (x_1,...x_i1+1, ...   x_i2+1, .... x_in+1,.... x_D) -> y + sum_k=1^d L^(ik-1)

%to do the above calculaten, first create a array with all directions
directions = nchoosek(1:D,d);
%These directions have to be translated to into additions:
% for example, if D=4 and d=2
% 1 2 -> [0 , 1,   L,   1+L]
% 1 3 -> [0 , 1,   L^2, 1+L^2]
% 1 4
% 2 3
% 2 4
% 3 4 -> [0 , L^2, L^3, L^2+L3]

% first, create a binary matrix of size 2^d, from which we can create the
% 2^d different vertices
M = [];
for i=0:d-1
    M = [kron([0;1],ones(2^i,1)),repmat(M,[2,1])];
end
%for example, if d=2;
% 0 0
% 0 1
% 1 0
% 1 1
% It is used to decide what directions to take

if isempty(M)
    %the case when d=0, ie. 0-cells
    shifts = 0;
else
    shifts = L.^(directions-1)*M';
end
%[1, L                            [0, L  ,1 ,  1+L
% 1, L^2    x  [ 0 0 1 1      =    0, L^2,1 ,  1+L^2
% 1, L^3         0 1 0 1 ]
%....                              ....
% L^2, L^3]                        0,L^2,L^3,L^2+L^3]

% A matrix indicating which units cells we should consider
% we will only consider unit cells starting at some vertix with subscript
% [x1, x2....xD] with xi < L for all i.
M = padarray(ones([(L-1)*ones(1,D-1),f*(L-2)+1]),ones(1,D),0,'post');
M = M(:);

k = nchoosek(D,d);
%number of cells per unit cell;

%initialization, we have a grid of unit cells of length 
%  L-1          in spatial directions
%  (f*(L-2)+1)  in time directions
% and k cells per unit cell
dCells = zeros((f*(L-2)+1)*(L-1)^(D-1)*k,2^d);


r=1;
for i=1:length(M)
    % M checks if we want to take this unit cell.
    if M(i)==1
       dCells(k*(r-1)+1:k*r,:) = i + shifts;        
       r=r+1;
    end
    
end

end