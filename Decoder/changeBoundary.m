function Cells = changeBoundary(Cells,L,D,direction,boundarytype,RoughSmooth,f)
%either makes 'pre' boundaries rough or 'post' boundaries smooth

% the 'post' boundaries are those located at some coordinate x_j = L
% the 'pre'  boundaries are those located at some coordinate x_j = 1

% Cells         for a description see 'createCells.m'
% L             size of system
% D             dimension (including time)
% direction     determines the coordinate j
% boundarytype  is either 'post' or 'pre'
% RoughSmooth   is either 'R' or 'S' (for rough or smooth)
% f             a factor indicating if time direction is larger than
%               spatial directrions

% note that the RoughSmooth variable is inprinciple redundant because we
% only use it for the cases (boundarytype,RoughSmooth) = 
% 'pre' - 'R'
% 'post' - 'S'

% create a binary direction vector d of size D whose only non-zero entry is d_j
d = zeros(1,D);
d(direction)=1;

% V is a [L,L, ... L(f)] dimensional (binary) array of rank D, where
% V(x1,x2,.... xD) is only 1 if the vertix with subscript (x1,x2,... xD) is
% part of the boundary which needs to be changed
V = zeros([L*ones(1,D-1),f*(L-2)+2]-d);
V = padarray(V,d,1,boundarytype);

% find the linear indices of the vertices of the boundary
% linear indices are related to subscripts as i = 1 + \sum (x_i-1) L^(i -1)
V = find(V(:)>.5);

%number of cells to be checked
nCells = size(Cells,1);
% a boolean array, storing which cells will be removed.
keepCells = true(nCells,1);

% for each cell, check if it needs to be removed
for i=1:nCells
    switch RoughSmooth
        case 'R'
            % for rough boundary, any cell is removed if it contained in
            % the boundary
            removeCell = all(ismember(Cells(i,:),V));  
        case 'S'
            % for smooth boundary, any cell is removed if it has overlap
            % with the boundary
            removeCell = any(ismember(Cells(i,:),V));
    end
    keepCells(i) = ~removeCell;
end

%keep only those cells which are not removed.
Cells = Cells(keepCells,:);

end 
    
    
    
    