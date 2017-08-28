function [fq,weight] = findSurface(stabilizers,qubits,syndrome,SQ,wq,id)
global opt_count Tune

%optimize per pbox

%INPUTS
% stabilizers   boolean vector indicating which stabilizers are part of
%               the box
% qubits        boolean vector indicating which qubits are part of
%               the box
% syndrome      full syndrome info
% SQ            parity check matrix
% wq            log probability of flipping a qubit
% id            unique identifier for each box

% total number of qubits and stabilizers in a box
nq = sum(qubits);       
ns = sum(stabilizers);

% syndrome info and parity check matrix of a box
syndrome_box = syndrome(stabilizers);
SQ_box = SQ(stabilizers,qubits);

% gurobi minimizes an objective, our decoding problem aims to maximize a
% certain probability. This is adjusted with a minus sign.
wq_box = -wq(qubits);

% only run integer linear program if there are syndromes contained in the
% box
if sum(syndrome_box)>.5 || strcmp(Tune,'tune');

    %weight vector, the last "ns"-entries are zero since these are
    %dummy variables simply used to translate the mod 2 problem to an
    %integer problem
    f = [wq_box;zeros(ns,1)];

    %number of integer constraints
    intcon = 1:(nq+ns);

    %constrain matrix
    Aeq = [SQ_box  2*eye(ns) ];       
    beq = syndrome_box;

    %bounds:
    %bounds on qubits are simply between 0 and 1
    %bounds on dummy variables are betweeen -4 and 0 since their are
    %max 8 qubits per stabilizer (in 5 dimenions, in general their are
    %2(D-1) qubits per stabilizer )
    lb = [zeros(nq,1);-4*ones(ns,1)];
    ub = [ones(nq,1);zeros(ns,1)];

    %optimizing
    tic; 
    [f_out,weight,ef] = intlinprogGUROBI(f,intcon,[],[],Aeq,beq,lb,ub,id);
    t = toc;

    %print out a warning if optimization failed
    if ef~=1
        fprintf('!!!!!!!!!!!!!!____%d\t%f\n' ,ef,min(f))
    end

    %print out optimization run time when tuning the model
    if strcmp(Tune,'tune')
        opt_count = opt_count+1;
        fprintf('%d\t%d\t%d\t%d\t%f\n',opt_count,size(Aeq,1),size(Aeq,2),sum(beq),t)
    end
else
    % trivial solution if no syndrome is in the box
    f_out = zeros(nq+ns,1);
    weight = 0;
end

%only output first nq variables (not the dummy variables)
fq = f_out(1:nq);    
end