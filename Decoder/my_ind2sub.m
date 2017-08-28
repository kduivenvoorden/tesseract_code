function coordinates = my_ind2sub(L,indices,D,f)
% the matlab function ind2sub, adjusted to my own needs

% it outputs a single matrix, instead of D matrices one for each coordinate
% indices, an array (of size [L1,... Ln] of indices
% coordinates, an array of size [L1,... Ln, D] of the subscripts of the
% indices

% the indices are assumed to be some coordinate of a [L, L, ... L, fL]
% sized hypercube

% the dimension can be deduced from the global variable rs-global, when not
% given
global rs_global

% default value of f
if nargin == 3
    f = 1;
end

% default value of D and f
if nargin == 2
    D = length(rs_global);
    f = 1;
end

% default value of D
if isempty(D)
    D = length(rs_global);
end


% dimensions of the input matrix indices
ind_size = size(indices);
% reshape the input matrix into a vector, to use the ind2sub function of
% matlab
indices = indices(:);

%depending on the dimension D, use the appropriate output syntax for the
%matlab function ind2sub
switch D
    case 2
        [x,y] = ind2sub([L,f*L],indices);
        coordinates = [x,y];
    case 3
        [x,y,z] = ind2sub([L,L,f*L],indices);
        coordinates = [x,y,z];
    case 4
        [x,y,z,w] = ind2sub([L,L,L,f*L],indices);
        coordinates = [x,y,z,w];
    case 5
        [x,y,z,w,v] = ind2sub([L,L,L,L,f*L],indices);
        coordinates = [x,y,z,w,v];
        
end
%reshape the subscripts into the same form as the input matrix
coordinates = reshape(coordinates,[ind_size,D]);
end