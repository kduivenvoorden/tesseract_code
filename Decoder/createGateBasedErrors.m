function [syndrome, readout_Q_end,fq] = createGateBasedErrors(N,f,p,CNOTS_X,CNOTS_Z)

% create errors basesd on the circuit specified by a bunch of CNOTS and 

%INPUT
%N          size of lattice L=2^N+2
%f          possible extension of time direction (IS NOT IMPEMENTED IN THIS
            %FUNCTION)
%p          failure parameter
%CNOTS_X    cell specifying all CNOT gates for X measurement for each round
%CNOTS_Z    cell specifying all CNOT gates for Z measurement for each round

%OUTPUT
%syndrome       value of all erroneous stabilizer measurements
%readout_Q_end  true flipped qubit at end of circuit, used to check of
                %correction was ok
%fq             check which errors occured in a (space,time) lattice
                %used for debugging purpose and to create an efficte failure prob per qubit/ stab
                %which can be used to optimize the decoder
                %syndrome should (more or less) be the boundary of fq

%Global parameters
global rs_global lattice_st lattice_s
%rs_global      rough and smooth boundaries
%lattice_st     lattice in space time
%lattice_s      lattice, no time direciont

%% initialize

%boundaries of lattice, disregarding time
rs = rs_global(1:end-1);
%size of lattice
L=2^N+2;
%dimension of lattice
D=length(rs);
%definition of edges, X-stab
SX_def = lattice_s.dCells_all{2,N+1};
%definition of faces , qubits
qubit_def = lattice_s.dCells_all{3,N+1};
%definition of cubes , Z-stab
if D > 2
    SZ_def = lattice_s.dCells_all{4,N+1};
else
    SZ_def = [];
end;

% different probabilites
%ancilla preperation
p_prep  =p;
%ancilla measrement
p_meas  =p;
%qubit idiling , depolarizing 
p_idle_Q=p;
%X-stab ancilla idiling , depolarizing 
p_idle_X=p;
%Z-stab ancilla idiling , depolarizing 
p_idle_Z=p;
%CNOT gate error
p_cnot  =p;

%number of qubits
nq = size(qubit_def,1);
%number of X-stabilizers
nsX = size(SX_def,1);
%number of Z-stabilizers
nsZ = size(SZ_def,1);




%matrix, keeping track of errors on qubits, either pauli X or pauli Z
e_Q = zeros(nq,2);
%first collumn: X error happend
%second collumn: Z error happend

%matrix, keeping track of all outcomes of X-stab measurement of L-1 rounds
readout_X = zeros(nsX,L-1);
%matrix, keeping track of all outcomes of Z-stab measurement of L-1 rounds
readout_Z = zeros(nsZ,L-1);
%matrix, keeping track of all errors on qubits, either pauli X or pauli Z,
%also knowing during which round the error happend, this is used to
%calculate fq
qubitflip_all = zeros(nq,2,L-1);

%% appliciation of circuit errors

%in total there are L-1 rounds, L-2 noisy rounds and last round is perfect
for time_round = 1:L-2
    %ancilla initiation, 
    %ancilla qubits have no errors to start of with
    e_SX=zeros(nsX,2);
    e_SZ=zeros(nsZ,2);
    %depeding on the type of ancilla, an extra X or Z pauli error occurs
    %with probability X
    %|0> state for Z-stabilizers -> X pauli errors can occur
    %|+> state for X-stabilizers -> Z pauli errors can occur
    e_SX(:,2) = rand(nsX,1)<p_prep;
    e_SZ(:,1) = rand(nsZ,1)<p_prep;
    
    %during ancilla creation, other qubits idle 
    e_Q = idling(e_Q,p_idle_Q);
      
    %store qubit errors
    qubitflip_all(:,:,time_round) = e_Q;
    
    %Application of CNOT gates
    for i=1:2*D
        %Perfect application of CNOT gates is simply defined by propagating
        % pauli Z from control to target and pauli X from target to control
        
        % X-stab
        %qubits are target and obtain extra X error if corresonpind X-stab
        %ancilla does
        e_Q(CNOTS_X{i}(:,1),1) = e_Q(CNOTS_X{i}(:,1),1) + e_SX(CNOTS_X{i}(:,2),1);
        %ancilla are control and obtain extra Z error if corresonpind qubit does
        e_SX(CNOTS_X{i}(:,2),2) = e_SX(CNOTS_X{i}(:,2),2) + e_Q(CNOTS_X{i}(:,1),2);
        
        % Z-stab
        %qubits are control and obtain extra Z error if corresonpind Z-stab
        %ancilla does
        e_Q(CNOTS_Z{i}(:,1),2) = e_Q(CNOTS_Z{i}(:,1),2) + e_SZ(CNOTS_Z{i}(:,2),2);
        %ancilla are target and obtain extra X error if corresonpind qubit does
        e_SZ(CNOTS_Z{i}(:,2),1) = e_SZ(CNOTS_Z{i}(:,2),1) + e_Q(CNOTS_Z{i}(:,1),1);
        
        % error after CNOT of X-stabilizer
        % qubits on which CNOTS acts
        ind_q = CNOTS_X{i}(:,1);
        % stabilizers on which CNOTS acts
        ind_s = CNOTS_X{i}(:,2);
        %apply error
        [e_Q(ind_q,:),e_SX(ind_s,:)] = CNOT_error(e_Q(ind_q,:),e_SX(ind_s,:),p_cnot);
        
        % error after CNOT of Z-stabilizer
        % qubits on which CNOTS acts
        ind_q = CNOTS_Z{i}(:,1);
        % stabilizers on which CNOTS acts
        ind_s = CNOTS_Z{i}(:,2);
        %apply error
        [e_Q(ind_q,:),e_SZ(ind_s,:)] = CNOT_error(e_Q(ind_q,:),e_SZ(ind_s,:),p_cnot);
        
        %on all qubits on which no CNOTS act, apply idling gates 
        %list of all qubits
        ind =  1:nq;
        %remove those from list on which CNOTS acted
        ind([CNOTS_X{i}(:,1);CNOTS_Z{i}(:,1)]) = [];
        %apply depolarizing noise
        e_Q(ind,:) = idling(e_Q(ind,:),p_idle_Q);
        
        %list of all X-stab ancillas
        ind = 1:nsX;
        %remove those from list on which CNOTS acted
        ind(CNOTS_X{i}(:,2)) = [];
        %apply depolarizing noise
        e_SX(ind,:) = idling(e_SX(ind,:),p_idle_X);
        
        %list of all Z-stab ancillas
        ind = 1:nsZ;
        %remove those from list on which CNOTS acted
        ind(CNOTS_Z{i}(:,2)) = [];
        %apply depolarizing noise
        e_SZ(ind,:) = idling(e_SZ(ind,:),p_idle_Z);
    end
    
    %ancilla readout
    %with no errors occuring, the ancilla should simply be in the 
    %|+> state for X-stabilizers
    %|0> state for Z-stabilizers
    
    %first add extra noise to similate noisy ancilla measurement
    e_SX(:,2) = e_SX(:,2) + (rand(nsX,1)<p_meas);    
    e_SZ(:,1) = e_SZ(:,1) + (rand(nsZ,1)<p_meas);
    
    % for X-stab, odd number of pauli-Z error implies an excitation
    % for Z-stab, odd number of pauli-X error implies an excitation
    readout_X(:,time_round) = e_SX(:,2); 
    readout_Z(:,time_round) = e_SZ(:,1); 
    
    %idling time of qubits during ancilla readout
    e_Q = idling(e_Q,p_idle_Q);
    
    %last round (L-1) is perfect and is extra flippeds are introduced, 
    %nevertheless qubitflip_all(:,:,L-1) \neq qubitflip_all(:,:,L-2) due to
    %error creation and propagation within the L-2 round
    if time_round == L-2
        qubitflip_all(:,:,L-1) = e_Q;
    end
    
end

%since the final (L-1) readout is perfect one simply has to check the
%accumulation of pauli Z-errors on qubits, to figure out X-stab outcome
readout_X(:,end) = lattice_s.SQX_all{N+1}*e_Q(:,2); 
%and accumulation of pauli X-errors on qubits, to figure out Z-stab outcome
%(which is not used....)
readout_Z(:,end) = lattice_s.SQZ_all{N+1}'*e_Q(:,1); 
%final readout of accumulation of aal pauli Z errors to check the
%performance of the decoder
readout_Q_end = e_Q(:,2);

%% syndrome calculation
% stabilizer measurements of all rounds have to be translated to a syndrome
% of a space-time lattice. Moreover, syndrome corresponding to edges
% oriented in the time direction have to be put in place.

%space time lattice
%definition of edges, X-stabilizers
SX_def_all      = lattice_st.dCells_all{2,N+1};
%definition of faces, qubits
qubit_def_all   = lattice_st.dCells_all{3,N+1};
% stabilizers incident to a vertix
VSX             = lattice_st.VSX_all{N+1};
% parity check matrix of X-stab
SQX             = lattice_st.SQX_all{N+1};

% initialize syndrome
syndrome = zeros(size(SX_def_all,1),1);

%syndromes are on edges oriented either in time or in space direction,
%those oriented in space direction can simply be obtained by the difference
% of measurement at time t and t-1. These don't form closed loops due to
% syndrome in time direction. These will be fixed by closing the loop.

% indices of syndromes which will not be changed to close the loop of
% syndromes
fixed_ind = false(size(SX_def_all,1),1);




% initialize fliped-qubits
fq = zeros(size(qubit_def_all,1),1);

%fq are on face oriented either in space-time or in space direction,
%those oriented in space direction can simply be obtained by the difference
% of qubit flip at time t and t-1. Those oriented in in space-time
% direction will be determined by by checking the diffence between the boundary of qubit
% flip and syndrome (which is accounted to syndrome measurement failure)

% indices of faces which will not be flipped to account for measurement failure
fixed_ind_qub = false(size(qubit_def_all,1),1);

% for each measurement round
for time_round = 1:L-1
    %edges, syndrome in space direction in a time slice
    hor_ind_syn = all(ismember(SX_def_all,SX_def+(time_round-1)*L^D),2);
    %faces, qubits in space direction in a time slice
    hor_ind_qub = all(ismember(qubit_def_all,qubit_def+(time_round-1)*L^D),2);
    
    %calculation of fq and syndrome
    if time_round ==1
        %for first measurement, syndrome(1) = readout(1)
        syndrome(hor_ind_syn) = readout_X(:,time_round);
        %for first measurement, fq(1) = qubitflip(1)
        fq(hor_ind_qub) = qubitflip_all(:,2,time_round);
    else
        %for all other measurement, it is the difference
        syndrome(hor_ind_syn) = readout_X(:,time_round)-readout_X(:,time_round-1);
        fq(hor_ind_qub) = qubitflip_all(:,2,time_round)-qubitflip_all(:,2,time_round-1);
    end
    %edges and faces in space direction are fixed and will not be altered
    %to close the erroneous syndrome loop
    fixed_ind(hor_ind_syn) = true;
    fixed_ind_qub(hor_ind_qub) = true;
    
    %vertices in the time-slices
    vertices = any(VSX(:,hor_ind_syn),2);
    

    % syndrome joing time slices
    if time_round<L-1
        % find the edges in the time direction by checking taking the
        % vertices in a time slice, considering the incident edges, and not
        % considering those edges which are in space direction (or previous
        % time slice). 
        % ver_ind does not contain repeated indices since
        % there is only one edge per vertix sticking in positive time
        [~,ver_ind,~] = find(VSX(vertices,:) & repmat(~fixed_ind',[sum(vertices),1]));
        % calculate VSX*syndrome, taking only the fixed syndrome, ie. those
        % in space. This should be equal to the corresponding syndrome in
        % time direction in order for the loop to close
        syndrome(ver_ind) = VSX(vertices,fixed_ind)*syndrome(fixed_ind);
        
        %consider all faces incident to edges in a space direction in a
        %time slice
        [~,ver_ind_qub,~] = find(SQX(hor_ind_syn,:)) ;
        %remove all faces wchi are all ready fixed, ie. those in space
        %direction or those in previous time direction
        ver_ind_qub(ismember(ver_ind_qub,find(fixed_ind_qub)))=[];
        %what remains are only those faces in space-time direction
   
        %flip those qubits if the space syndrome cannot be accounted
        %for by space flips
        fq(ver_ind_qub) = mod(syndrome(hor_ind_syn) +SQX(hor_ind_syn,fixed_ind_qub)*fq(fixed_ind_qub),2);
    end
    
    % fix all edges and faces in the time direction
    fixed_ind(ver_ind) = true;
    fixed_ind_qub(ver_ind_qub) = true;
end

%syndrome is either 0 or 1
syndrome = mod(syndrome,2);


end


function [X,Y] = CNOT_error(X,Y,p)
% X     N by 2  matrix containing (X,Z) errors of first  qubit
% Y     N by 2  matrix containing (X,Z) errors of second qubit
% p     failure probablity
    if p>0
        %matrix containing the different options
        SQ =[...
     0     0     0     0;...
     0     0     0     1;...
     0     0     1     0;...
     0     0     1     1;...
     0     1     0     0;...
     0     1     0     1;...
     0     1     1     0;...
     0     1     1     1;...
     1     0     0     0;...
     1     0     0     1;...
     1     0     1     0;...
     1     0     1     1;...
     1     1     0     0;...
     1     1     0     1;...
     1     1     1     0;...
     1     1     1     1];
        x =rand(size(X,1),1);


        n = 15;
        x(x<1-p) = (1-p)-p/(2*n);
        %index, labeling the option 
        y = 1+round((x-(1-p))/p*n+.5);
        %add error to the already exsiting noise
        X = X+SQ(y,1:2);
        Y = Y+SQ(y,3:4);
    end
end


function X = idling(X,p)
% X     N by 2  matrix containing (X,Z) errors
% p     failure probablity

% depolarizing noise,
% prop p/3, X-error (1,0)
% prop p/3, Y-error (1,1)
% prop p/3, Z-error (0,1)
% prop 1-p, I-error (0,0)
    if p>0
        %matrix containing the different options
        SQ = [0,0;0,1;1,0;1,1];
        x =rand(size(X,1),1);
        %index, labeling the option 
        y = 1+(x<p/3) + (x<2*p/3) + (x<3*p/3);
        %add error to the already exsiting noise
        X = X + SQ(y,:);
    end
end