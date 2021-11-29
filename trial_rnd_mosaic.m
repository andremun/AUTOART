tid = str2double(getenv('SLURM_ARRAY_TASK_ID'));
disp(['Trial number: ' num2str(tid)]);

a = 1:1000:1e6;
b = a(2)-1:1000:1e6;
range = a(tid):b(tid);

datadir = './autoart_1e6_cost/';
getfromfile = @(filename,varname) getfield(load(filename,varname),varname);

idx = getfromfile([datadir 'img_idx_1e6.mat'],'idx');
ntries = size(idx,2);
J = NaN.*ones(5,ntries);
J(4,:) = getfromfile([datadir 'result_gen_rand_mosaics_E0.mat'],'J');

Jaux = test_random_mosaics_clust(datadir,idx(:,range),J(:,range));
save([datadir 'result_gen_rand_mosaics_TID' num2str(tid) '.mat'],'Jaux');