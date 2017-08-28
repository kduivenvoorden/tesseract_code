function peform_runs(ntrials,Ls,ps,f,file_id,m_or_s,type)


%peform multiple runs for a certain parameter set
%INPUT
%ntrials    number of runs
%Ls         sytem sizes (could be vector)
%ps         probabilities (could be vector)
%f          time extension
%file_id    a string used for saving data
%m_or_s     a string used decide between multi or single thread
%type       a string used 'fem' or 'gb' or 'perf_meas' or '3d' determining the error model 
%           'fem' = femenological, independent X/Z
%           'gb'  = gated based


%make sure matlab has acces to gurobi and all other .m files
current_path = pwd;
cd /usr/local/gurobi605/linux64/matlab
gurobi_setup
cd(current_path)
addpath(strcat(pwd,'/Decoder'))
%code to set the LICENSE FILE environment variable (of the shell matlab is running in) to the right directory
%system('setenv GRB_LICENSE_FILE /mount/vol1/scratch/work/kasperd/gurobi.lic');

%check if the directories exist
dirs = {'Cells','Maps','CNOT',strcat('results_local_',type)};
for i=1:length(dirs)
    if exist(dirs{i},'dir')==0
        mkdir(dirs{i})
    end
end

%default values
if nargin ==0
    ntrials = 1;
    Ls = 4;
    ps = 0.04;
    f=1;
end



% length of parameter set
np = length(ps);
nL = length(Ls);

% For all system sizes
for iL = 1:nL
    L=Ls(iL);
    N = log2(L-2);
    %for all probabilities
    for ip=1:np
        p = ps(ip);

        %if multithread ...
        switch m_or_s
            case 'm'
                % ... create a pool for parallel threads
                my_pool = gcp();
                x = my_pool.NumWorkers;
                fprintf('number of workers: %d\n',x)
                %number of trials per thread
                trial_per_pools = ceil(ntrials/x);
                % for each thread
                parfor pool_id=1:x
                    run(p,N,f,trial_per_pools,file_id,pool_id,type)
                end % for thread
            case 's'
                    run(p,N,f,ntrials,file_id,1,type)
        end
        
    end % for probability
end % for system size
    
end    

function run(p,N,f,trial_per_pools,file_id,pool_id,type)
global rs_global visualizes Tune lattice_st lattice_s

%fixed parameters
% rough/smooth boundaries 0 smooth, 1 rough
% two rough boundary needed in order for a surface like X-logical
rs          = [0;0;1;1];



%Gurobi optimization
Tune        = 'no';
%output visualization of the decoder performance
visualizes = false;

%extra variable, keeping track if measurements are perfect
switch type
    case 'fem'
        type_struct.noise_model = 'fem';
        type_struct.meas_noise  = true;
    case 'perf_meas'
        type_struct.noise_model = 'fem';
        type_struct.meas_noise  = false;
    case 'gb'
        type_struct.noise_model = 'gb';
        type_struct.meas_noise  = true;
    case '3d'
        type_struct.noise_model = 'fem';
        type_struct.meas_noise  = false;
        rs = [0;1;1];
end






%initialze the lattice
lattice_s = createMapsCells(N,rs,f);

%depending if measurement is perfect
if type_struct.meas_noise
    %add time dimension for space-time lattice
    rs = [rs;0];
    lattice_st = createMapsCells(N,rs,f);   
else
    %space-time lattice is same as space-lattice
    lattice_st = lattice_s;
end

%global variable storing the type of boundaries, used by the decoder
rs_global = rs;    

%for all trials
for trial = 1:trial_per_pools
    fprintf('%d ',trial)

    try
        %decode
        result = localDecoder(p,N,f,type_struct); 

        %save the result
        fid = fopen(strcat('results_local_',type,'/results_',file_id,sprintf('_pp%d_D%d_f%d',pool_id,length(rs),f)),'a');
        fprintf(fid,datestr(datetime));
        fprintf(fid,'\t L=%02d p=%.5f t=%01d\n',2^N+2,p,result);
        fclose(fid);
    catch em
        em
    end
end

end
