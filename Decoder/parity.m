function p = parity(vertices,L)
    % determine the number of even subscripts s_i of the index x 
    % x  = (s_i-1) L^(i-1) +1 
    
    % vertices can be an D-dim array, p will be an array of the same size
    
    % determine subscripts of the indices
    p = my_ind2sub(L,vertices);
    
    %the origin is (1,1,1....,1), the function is created such that this
    %has parity 0
    p=p-1;
    p=mod(p,2);
    
    %sum over the last dimenions, this dimension stores the different
    %subscript of a certain index
    p = sum(p,length(size(p)));
    
end