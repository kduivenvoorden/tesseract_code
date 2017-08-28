%some examples of the workings of the decoder

N=      1;      %number of runs 
L=      6;      %system size (choose between 3,4,6,10)
p=      0.07;   %probability (threshold is approx 7.3%)
f=      1;      %ratio between time and space direction, is always set to 1
%5th argument is a string used for data-saving
cores = 's';    %single machine, set to 'm' for parrallel programming

%% perfect measurement fenomenological noise
peform_runs(N,L,p,f,'test',cores,'perf_meas')
cores = 'm'; %Always run on single machine before proceding to parrallel 
             %to avoid conflicts when creating adjacency matrices
N=10;
peform_runs(N,L,p,f,'test',cores,'perf_meas')

%% 3d with fenomenological noise
peform_runs(1,L,p,f,'test','s','3d')
peform_runs(N,L,p,f,'test','m','3d')

%% 4d with noisy measurement
L=      4;
p=      0.04;
peform_runs(1,L,p,f,'test','s','fem')
peform_runs(N,L,p,f,'test','m','fem')

%% 4d with gate based measurement
p=      0.003;
peform_runs(1,L,p,f,'test','s','gb')
peform_runs(N,L,p,f,'test','m','gb')

%% make pictues
localDecoder