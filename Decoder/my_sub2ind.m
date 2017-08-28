function indices = my_sub2ind(L,coordinates,D,f)
% the matlab function sub2ind, adjusted to my own needs

% it has as input single matrix, instead of D matrices one for each coordinate
% coordinates, an array of size [L1,... Ln, D] of the subscripts of the
% indices
% indices, an array (of size [L1,... Ln] of indices

% the indices are assumed to be some coordinate of a [L, L, ... L, fL]
% sized hypercube

% the dimension can be deduced from the global variable rs-global, when not
% given
global rs_global

% default value of f
if nargin == 3
    f = 1;
end

% default value of f and D
if nargin == 2
    D = length(rs_global);
    f = 1;
end

% default value of D
if isempty(D)
    D = length(rs_global);
end

% dimensions of the input matrix coordinates
coordinate_size = size(coordinates);
% reshape the input 2-dim array, to use the sub2ind function of
% matlab
coordinates = reshape(coordinates,[],D);

%depending on the dimension D, use the appropriate input syntax for the
%matlab function sub2ind
switch D
    case 2
        indices = sub2ind([L,f*L],coordinates(:,1),coordinates(:,2)); 
    case 3
        indices = sub2ind([L,L,f*L],coordinates(:,1),coordinates(:,2),coordinates(:,3));
        
    case 4
        indices = sub2ind([L,L,L,f*L],coordinates(:,1),coordinates(:,2),coordinates(:,3),coordinates(:,4));
    case 5
        indices = sub2ind([L,L,L,L,f*L],coordinates(:,1),coordinates(:,2),coordinates(:,3),coordinates(:,4),coordinates(:,5));
end
%reshape the indices into the same form as the input matrix
indices = reshape(indices,coordinate_size(1:end-1));
end