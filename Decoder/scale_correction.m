function sc  = scale_correction(correction,N,f)
% scales up qubit corrections, from a smaller to a larger lattice

% corrections   a binary vector of length nQ
% N             size of current lattice L = 2^N+2
% f             possible extension in time dimension
global lattice_st


L=2^N+2;


%Qubits are faces (2 cells)
qubit_def = lattice_st.dCells_all{3,N+1};

%vertices of qubits that need to be corrected
correction_def = qubit_def(mod(correction,2)==1,:);

%subscripts of the vertices of qubits that need to  be corrected
coordinates = my_ind2sub(L,correction_def,[],f);

%scale up subscripts, to match a 2 by 2 face of the larger lattice
coordinates = (coordinates-1)*2+1;

%take the average of the coordinates of the 4 vertices, giving the coordinate of the
%center of the 2 by 2 face
coordinates = mean(coordinates,2);

%size of larger lattice
new_L = 2^(N+1)+2;

%the index of the center of a 2 by 2 face, of 4 qubits that need to be
%corrected
vertex_def = my_sub2ind(new_L,coordinates,[],f);

%Qubits of larger lattice 
new_qubit_def = lattice_st.dCells_all{3,N+2};

%the qubits that need to be corrected are those that have a vertex being
%one of the centers of a 2 by 2 face and more over overlap with the coarse
%grained lattice. The latter condition can be checked by means of the
%parity of the vertices. A face being on the coarse grained lattice has
%vertices with parity [0,1,1,2]

%qubits bigger lattice which are also on coarse grained lattice
qubit_on_course_grained = sum(parity(new_qubit_def,new_L),2)==4;

%qubits bigger having a vertix being the center of a 2 by 2 face
qubits_on_to_be_corrected_22face  = any(ismember(new_qubit_def,vertex_def),2);


sc = qubits_on_to_be_corrected_22face & qubit_on_course_grained;

end