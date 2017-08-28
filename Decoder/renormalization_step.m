function [new_syndrome,new_wq,correction] = renormalization_step(syndrome,wq,N,f)
%
%moves syndrome to coarse grained lattice

%INPUTS
% syndrome  binary vector of length ns (number of stabilizers)
% wq        vector of length nq (number of qubits), wq_i is the probabilty
%           that the qubit i flipped
%           q_i = log(pf/pn) where pf is probabilty of flip, and pn is
%           probability of no flip
% N         size of lattice L = 2^N+2
% f         possible extension of time direction

% OUTPUTS
% new_syndrome  binary vector of length ns' (number of stabilizers in
%               coarse grained lattice
% new_wq        vector of length nq, (rescalled) probalities of qubit
%               corrections in coarse grained lattice
% correction    binary vector of length nq


 
%% Preperation
global rs_global visualizes lattice_st

% Dimension of lattice
D               = length(rs_global);
% stabilizer definitions
SQX             = lattice_st.SQX_all{N+1};
% stabilizers are edges
stabilizer_def  = lattice_st.dCells_all{2,N+1};
%qubits are faces
qubit_def       = lattice_st.dCells_all{3,N+1};
%highest dimenional cube
cube_def        = lattice_st.dCells_all{D+1,N+1};
% all vertices of the lattice which are relavent
vertices_def    = unique(qubit_def);
% size of lattice
L               = 2^N+2;
%number of qubits and stabilizers
[ns,nq]         = size(SQX);
%number of vertices
nv = size(vertices_def,1);
%correction applied to qubits
correction      = zeros(nq,1);

if visualizes
    %plot the lattice and the syndrome info
    visualize1(L,stabilizer_def,syndrome,f,D)
end
  
%% end step: solve a global interger program
if N==0    
    
    % used to save a certain program, for gurobi optimization 
    id =-1;
    [correction,~ ] = findSurface(true(ns,1),true(nq,1),syndrome,SQX,wq,id);
    
    new_syndrome = 0;
    new_wq = 0;
else
%% renormalization step: initialization
    
    %Parity is an important quantity to check the location of certain cells
    %relative to the course grained lattice. It is for example used to
    %determine which edges are also edges of the course grained lattice or
    %to determine the center of optimization boxes.
    vertices_parity = parity(vertices_def,L);
    %parity of vertices is either 0,1,...D
    %0 : also vertix of course grained (all even)
    %1 : on edge of course grained (one odd)
    %2 : on face of course grained (two odd)
    %3 : in cube of course grained (three odd)
    %...

    %devide edges into number of odd components
    stabilizer_parity = (sum(parity(stabilizer_def,L),2)-1)/2;
    %parity of stabilizer is either 0,1,...D-1
    %0 : also on edge of course grained 
    %1 : on face of course grained 
    %2 : in cube of course grained 
    %...

    qubit_parity = (sum(parity(qubit_def,L),2)-4)/4;
    %parity of qubit is either 0,1,...D-2
    %0 : on face of course grained 
    %1 : on cube of course grained 
    %...
    

    % number of boxes per stabilizer,
    % The lattice wil be devided into boxes, which could be overlapping. 
    % number_of_boxes_per_stabilizer is a vector of length ns keeping track of
    % how many boxes a stabilizer is contained in. This is important since if
    % all boxes containing a certain stabilizer are handled, that stabilizer must be
    % free of syndrome (unless it is also overlapping with the course
    % grained lattice).
    number_of_boxes_per_stabilizer = zeros(ns,1);
    for ivertix= 1:nv
       %boxes are identified by their center vertices, which have parity D
       if vertices_parity(ivertix)==D 
           %corresponding index of vertix
           ind_center = vertices_def(ivertix);
           %box contains all vertices which are distance 1 away from the
           %center vertix (in l_infty norm)
           box = cube_def(any(cube_def==ind_center,2),:);       
           
           %stabilizers contained in a certain box
           stabilizers_in_box = all(ismember(stabilizer_def,box),2)&stabilizer_parity~=0;
           %increase their number_of_boxes_per_stabilizer count
           number_of_boxes_per_stabilizer(stabilizers_in_box)  = number_of_boxes_per_stabilizer(stabilizers_in_box)+1;
       end
    end
  
%% renormalization step: optimization per box       

    %booleans to remove certain qubits from an integer program if they
    %are bounded by a stabilizer which will no longer be contained in any
    %to be optimized box, and as such, should stay free of syndrome
    stabilizer_count = zeros(ns,1);
    qubits_not_to_be_changed = false(nq,1);
    stabilizer_not_to_be_changed = false(ns,1);
    
    % new_syndrome will be the syndrome after correction (=syndrome
    % + SQX*correction)
    new_syndrome = syndrome;
    
    % used for debugging puprposes
    box_count=0;
    
    for ivertix= 1:nv
        if vertices_parity(ivertix)==D
            box_count=box_count+1;

            % identifier for a certain box, used to save a certain program, for gurobi optimization 
            id = 10*ivertix+N;

            %box contains all vertices which are distance 1 away from the
            %center vertix (in l_infty norm)
            ind_center = vertices_def(ivertix);
            box = cube_def(any(cube_def==ind_center,2),:);

            %stabilizers in the integer program
            stabilizers_in_box = all(ismember(stabilizer_def,box),2)&stabilizer_parity~=0;
            %keep track of how many times a stabilizer is contained in a
            %to-be-optimized box
            stabilizer_count(stabilizers_in_box) = stabilizer_count(stabilizers_in_box)+1;
            % and flag them as soon as the last box containing a certain stabilizer is being optimized           
            stabilizer_not_to_be_changed(stabilizer_count>=number_of_boxes_per_stabilizer)=true;
            % this never holds for stabilizer, which are also part of the
            % coarse grained lattice. These need not be free of syndrome
            stabilizer_not_to_be_changed(stabilizer_parity==0) = false;

            %All qubits, overlapped by a certain stabilizer which is part of
            %the integer program, should also be part of the integer
            %program, unless it should not be flipped any more (since it
            %will create a syndrome not on the course grained lattice which
            %can no longer be removed).
            qubits_in_box = any(SQX(stabilizers_in_box,:))' & ~qubits_not_to_be_changed;       

            %perform the integer program
            [fq_box,~] = findSurface(stabilizers_in_box,qubits_in_box,new_syndrome,SQX,wq,id);    

            %update correction and the new_syndrome
            correction(qubits_in_box) = mod(correction(qubits_in_box)+fq_box,2);         
            new_syndrome = mod(syndrome+SQX*correction,2);

            %update probabilities, the probablity to flip it back is - the
            %original probability
            wq(qubits_in_box) = wq(qubits_in_box).*((-1).^fq_box);

            %dont change these qubits, since it will create a syndrome which 
            % can not be removed any more
            qubits_not_to_be_changed(any(SQX(stabilizer_not_to_be_changed,:) )) = true;
            
            if visualizes && N==2
                %visualize an example of single box
                visualize3(N,f,box_count,qubits_in_box,stabilizers_in_box)
            end
       end
    end
    if visualizes
        %visualize the applied correction and the change in syndrome
        visualize2(L,D, stabilizer_def,new_syndrome,syndrome,correction,f,stabilizer_parity,qubit_def)
    end
   
    %sparse version of correction (not sure if it is needed)
    [r,c,v] = find(correction);
    correction = sparse(r,c,v,size(correction,1),size(correction,2),nnz(correction));
     
%% renormalization step: recalculate weights of course grained qubits

    % A face of the course grained lattice consits of four (parallel,
    % overlapping) faces of the lattice. It is (uniquely) identified by the
    % midway vertex, ie a vertex with parity 2. All midway vertices
    % identify an face unless they are on the preboundary of rough edge.
    
    %subscript of vertices
    vertices = my_ind2sub(L,vertices_def,D,f);
    vertices = squeeze(vertices);
    
    %keep only those subscript belonging to a rough direction
    vertices_rough = vertices(:,rs_global==1);
    %if any of those is 1, that vertex will not identiy an edge
    vertices_rough_preboundary = any(vertices_rough==1,2);   
    
    %coarse grained qubits
    new_qubits = (~vertices_rough_preboundary) & (vertices_parity==2);
    %and their corresponding weights
    new_wq = zeros(nv,1);
    
    
    % qubits which are on part of a face of the coarse grained lattice
    qubits_on_face_all = qubit_def;
    qubits_on_face_all(qubit_parity~=0,:) = 0;
    
    for ivertex= find(new_qubits)';
        %index of the vertix identifying a coarse grained edge
        vertex_loc = vertices_def(ivertex);
        
        %The four qubits which make up the current coarse grained qubits
        qubits_2_newqubit = any(qubits_on_face_all==vertex_loc,2);
        
        %weight of the four qubits
        wq_part = wq(qubits_2_newqubit);
        
        % the weight of the coarse grained qubit is in principle the
        % sum of the weights of the the four qubits (flipping a coarse
        % grained qubit is equivalent to flipping these four qubits).
        %new_wq(ivertex) = ((-1).^correction(qubits_2_newqubit))'*wq_part;
        new_wq(ivertex) = sum(wq_part);
    end
    


%% renormalization step: course grain lattice
    
    %defintions of the stabilizers on the coarse grained lattice
    new_stabilizer_def = lattice_st.dCells_all{2,N};
    %defintions of the qubits on the coarse grained lattice
    new_qubit_def = lattice_st.dCells_all{3,N};
    %size of the coarse grained lattice
    new_L = 2^(N-1)+2; 

    % stabilizers of the coarse grained lattice consist of pairs of edges of
    % the original lattice. The pairs are created by looking at the
    % vertices of these edges, one of them has parity 0 the other has
    % parity 1. The parity-1-vertices of two edges forming a pair are equal 
    
    % the parity of the vertices belonging of the different stabilizers
    x = sum(mod(my_ind2sub(L,stabilizer_def,D,f)-1,2),3);
    % only consider those that are part of the coarse grained lattice
    x(stabilizer_parity~=0,:)=0;
    
    % determine the vertix of the edges lattice which has parity 1
    y = stabilizer_def;
    y(x==0)=0;
    y = sort(y,2);
    y = y(:,end);
    
    % define pairs of edges, whenever their parity-1-vertices are the same.
    [~,~,iz] = unique(y);
    %if iz(n) = iz(m) = k then the edges n and m together form the k-th pair
    %all iz(n) = 1 for those n for which y(n) = 0. The first "pair" is not
    %a pair but merely is the group of zero entries of y.
    [stab_in_pair,~,pair] = find(iz-1);
    % put this info in a matrix
    PS = sparse(pair,stab_in_pair,1,max(pair),ns);
    
    % The next step is to elongate those "pairs" on rough boudaries that
    % consits of only one edge. First the definition of the lattice will be made larger
    % from L to L+1, by changing the indexing of the vertices
    stabilizer_def_2 = my_sub2ind(L+1,my_ind2sub(L,stabilizer_def,D,f),D,f);
    % next determine those values in the vector "pair" which only have a single
    % occurence
    PS_single = PS;
    PS_single(sum(PS,2)~=1,:)=0;
    [single_pair,single_stab,~] = find(PS_single); 
    % and determine their definition (in terms of vertices)
    single_pair_def = stabilizer_def_2(single_stab,:);
    
    % double the single pairs. Take a single pair, and shift it along its
    % orrientation [x, x+dx ] -> [x+dx,x+2dx]
    single_pair_partner_def = single_pair_def+ repmat(single_pair_def(:,end)-single_pair_def(:,1),1,2);
    
    % add these parners to the existing edges
    stabilizer_def_2 = [stabilizer_def_2;single_pair_partner_def];
    %and update the stab_in_pair and pair
    stab_in_pair = [stab_in_pair; ns+(1:length(single_pair))']; 
    pair = [pair; single_pair];
    %The pairs are now elongated
    
    %take only those edges which form pairs
    pair_def = stabilizer_def_2(stab_in_pair,:);
    %sort the pairs such that stabilizers forming the same pair are
    %consequtive
    [~,o]=sort(pair);
    pair_def = pair_def(o,:);
    % join the edges: [x, x+dx ] + [x+dx,x+2dx] = [x+2dx]
    pair_def = reshape(pair_def',4,[])';
    pair_def = [min(pair_def,[],2) max(pair_def,[],2)];
    % determine the corresponding vertices on a smaller lattice (by
    % deviding the coordinates by 2
    pair_coor = my_ind2sub(L+1,pair_def,D,f);
    pair_coor = (pair_coor-1)/2+1;
    new_pair_def = my_sub2ind(new_L,pair_coor,D,f);
    
    %The newly obtained list of edges E1 should coincide with
    %new_stabilizer_def, the edges of the coarse grained lattice E2. However,
    %the ordering of these list might not be equal. Hence in, the following
    % we find the map m such that E2(i) = E1(m(i))
    [~, a] = sortrows(new_stabilizer_def);
    [~, b] = sortrows(new_pair_def);
    [~,c]  = sort(a);
    coarsegrain_mapping = b(c);
    
    % PS maps the edges to pairs...
    new_syndrome = (PS*new_syndrome)>.5;
    % ... which coincide with edges of the coarse grained lattice
    new_syndrome = new_syndrome(coarsegrain_mapping);
    % these maps are here used to determine the syndrome on the coarse
    % grained lattice


    %qubitweight calculation-------------------------------------------
   
    % new_qubits is a boolean vector of size nv (number of vertices)
    % indicating the center of 2by2 faces being also qubits of the coarse
    % grained lattice
    % new_wq is a vector of size nv corresponding to their weights
    % Here, using the vertex indices we will determine how to connect
    % these 2by2 faces to the qubits of the faces of the coarse grained
    % lattice
    
    % subscript of center of 2by2 face
    coordinate = my_ind2sub(L,vertices_def(new_qubits),D,f);
    % the subscript consist of 2 even coordinates and the rest is odd. The
    % two even coordinates correspond to the orientation of the face
    % from [o,... e, ... e,...o] we construct
    %  [o,... e-1, ... e-1,...o] and
    %  [o,... e+1, ... e+1,...o]
    % ie. the two corners of the 2by2 face
    coordinate = coordinate-1;
    corners = repmat(zeros(size(coordinate)),[1,2,1]);
    corners(:,1,:) = coordinate-mod(coordinate,2);
    corners(:,2,:) = coordinate+mod(coordinate,2);

    % indices of the corners in the coarse grained lattice
    corners = corners/2+1;
    corners = my_sub2ind(new_L,corners,D,f);
    
    % these coincide with the 1st and 4th vertix of the faces in the coarse
    % grained lattice. The following will get the ordering right.
    [~, a] = sortrows(corners);
    [~, b] = sortrows(new_qubit_def(:,[1,4]));
    [~,c]  = sort(b);
    mapping = a(c);
    
    %ensure that the coarse grained cubits get the appropriate
    %weight
    new_wq = new_wq(new_qubits);
    new_wq = new_wq(mapping);
    
    
end

    


end

%% some undocumented functions used to create pictures
function visualize1(L,stabilizer_def,syndrome,f,D)

if D ==4
    [x,y,z,u] = ind2sub([L,L,L,L],stabilizer_def(mod(syndrome,2)==1,:));
    [x3,y3,z3,u3] = ind2sub([L,L,L,L],stabilizer_def);
    for t = 1:L-1
        intializeVisualization(L,f);
        pt = u(:,1)==t & u(:,2)==t;
        plot3(x(pt,:)'-1,y(pt,:)-1',z(pt,:)'-1,'r','LineWidth',5);
        pt = u3(:,1)==t & u3(:,2)==t;
        plot3(x3(pt,:)'-1,y3(pt,:)'-1,z3(pt,:)'-1,'r','LineWidth',.1,'LineStyle',':');     
    end
else
    intializeVisualization(L,f);
    [x,y,z] = ind2sub([L,L,L],stabilizer_def(mod(syndrome,2)==1,:));
    plot3(x'-1,y'-1,z'-1,'r','LineWidth',5);
    [x,y,z] = ind2sub([L,L,L],stabilizer_def);
    plot3(x'-1,y'-1,z'-1,'r','LineWidth',.1,'LineStyle',':'); 
    saveas(gcf,sprintf('syndrome_%d.pdf',round(log2(L-2))))
end


end

function visualize2(L,D, stabilizer_def,new_syndrome,syndrome,correction,f,stabilizer_parity,qubit_def)
if D==4
    sns = mod(new_syndrome,2)==1 & mod(syndrome,2)==1;
    [x1,y1,z1,u1] = ind2sub([L,L,L,L],stabilizer_def(mod(new_syndrome,2)==1&~sns,:));
    [x2,y2,z2,u2] = ind2sub([L,L,L,L],stabilizer_def(mod(syndrome,2)==1&~sns,:));
    [x3,y3,z3,u3] = ind2sub([L,L,L,L],stabilizer_def(sns,:));
    [x4,y4,z4,u4] = ind2sub([L,L,L,L],stabilizer_def(stabilizer_parity==0,:));
    [x5,y5,z5,u5] = ind2sub([L,L,L,L],qubit_def(correction==1,:));

    for t = 1:L-1
        intializeVisualization(L,f);
        pt = u1(:,1) == t & u1(:,2) == t;
        plot3(x1(pt,:)'-1,y1(pt,:)'-1,z1(pt,:)'-1,'Color',[1,0,0],'LineWidth',1.5*(D-1));

        pt = u2(:,1) == t & u2(:,2) == t;
        plot3(x2(pt,:)'-1,y2(pt,:)'-1,z2(pt,:)'-1,'Color',[0,1,0],'LineWidth',1.5*(D-1));

        pt = u3(:,1) == t & u3(:,2) == t;
        plot3(x3(pt,:)'-1,y3(pt,:)'-1,z3(pt,:)'-1,'Color',[0,0,1],'LineWidth',1.5*(D-1));

        pt = u4(:,1) == t & u4(:,2) == t;
        plot3(x4(pt,:)'-1,y4(pt,:)'-1,z4(pt,:)'-1,'Color',[0,0,0],'LineWidth',.1,'LineStyle',':');

        pt = u5(:,1) == t & u5(:,4) == t+1;
        plot3(x5(pt,[2,3])'-1,y5(pt,[2,3])'-1,z5(pt,[2,3])'-1,'Color',[0,0,0],'LineWidth',1*(D-1));
        pt = u5(:,1) == t-1 & u5(:,4) == t;
        plot3(x5(pt,[2,3])'-1,y5(pt,[2,3])'-1,z5(pt,[2,3])'-1,'Color',[0,0,0],'LineWidth',1*(D-1));

        pt = u5(:,1) == t & u5(:,4) == t;
        cor = fill3(x5(pt,[1,2,4,3])'-1,y5(pt,[1,2,4,3])'-1,z5(pt,[1,2,4,3])'-1,[0,1,0],'EdgeColor','None');
        alpha(cor,.2);
    end
else
    intializeVisualization(L,f);
    sns = mod(new_syndrome,2)==1 & mod(syndrome,2)==1;
    [x,y,z] = ind2sub([L,L,L],stabilizer_def(mod(new_syndrome,2)==1&~sns,:));
    plot3(x'-1,y'-1,z'-1,'Color',[1,0,0],'LineWidth',1.5*(D-1));

    [x,y,z] = ind2sub([L,L,L],stabilizer_def(mod(syndrome,2)==1&~sns,:));
    plot3(x'-1,y'-1,z'-1,'Color',[0,1,0],'LineWidth',1.5*(D-1));

    [x,y,z] = ind2sub([L,L,L],stabilizer_def(sns,:));
    plot3(x'-1,y'-1,z'-1,'Color',[0,0,1],'LineWidth',1.5*(D-1));

    [x,y,z] = ind2sub([L,L,L],stabilizer_def(stabilizer_parity==0,:));
    plot3(x'-1,y'-1,z'-1,'Color',[0,0,0],'LineWidth',.1,'LineStyle',':');

    [x,y,z] = ind2sub([L,L,L],qubit_def(correction==1,:));
    cor = fill3(x(:,[1,2,4,3])'-1,y(:,[1,2,4,3])'-1,z(:,[1,2,4,3])'-1,[0,1,0]);
    alpha(cor,.2);
    saveas(gcf,sprintf('correction_%d.pdf',round(log2(L-2))))
end
end

function visualize3(N,f,nbox,qubits_in_box,stabilizers_in_box)
global lattice_st
L = 2^N+2;
%qubits are surfaces
qubit_def = lattice_st.dCells_all{3,N+1};
stabilizer_def = lattice_st.dCells_all{2,N+1};
stabilizer_parity = (sum(parity(stabilizer_def,L),2)-1)/2;

if nbox ==1
    intializeVisualization(L,f);
    [x,y,z] = ind2sub([L,L,L],stabilizer_def(stabilizer_parity==0,:));
    plot3(x'-1,y'-1,z'-1,'Color',[1,0,0],'LineWidth',2,'LineStyle','-');
    loc =[2,2,2;4,2,2;2,4,2;4,4,2;2,6,2;4,6,2;2,2,4;4,2,4;2,4,4]-1;
    for i = 1:9
        text(loc(i,1),loc(i,2),loc(i,3),sprintf('%d',i),'Color',[0,1,0]','FontSize',40)
    end
end

if nbox ==9
    [x,y,z] = ind2sub([L,L,L],stabilizer_def(stabilizers_in_box,:));
    plot3(x'-1,y'-1,z'-1,'k:','LineWidth',3);
    [x,y,z] = ind2sub([L,L,L],qubit_def(qubits_in_box,:));
    cor = fill3(x(:,[1,2,4,3])'-1,y(:,[1,2,4,3])'-1,z(:,[1,2,4,3])'-1,[0,0,0]);
    alpha(cor,.1);
    saveas(gcf,'box9.png')
end


end
