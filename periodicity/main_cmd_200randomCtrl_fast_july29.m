% matlab -nodesktop -r main_cmd_sep15_m0m1_rand_1sp.m -logfile run.log
% nohup matlab -nodesktop -r "try, main_cmd_sep15_m0m1_rand_1sp_v2_all_part2;end;quit" -logfile run1013.log
% nohup matlab -nodesktop -r "try, main_cmd_sep15_m0m1_rand_1sp_v3_dec15;end;quit" -logfile run1215.log

% matlab -nodesktop -r "main_cmd_200randomCtrl_fast_apr21"
% nohup matlab -nodesktop -r "try, main_cmd_200randomCtrl_fast_apr21;end;quit" 
% nohup matlab -nodesktop -r "try, main_cmd_200randomCtrl_fast_jun13;end;quit" 
% nohup matlab -nodesktop -r "try, main_cmd_200randomCtrl_fast_july29;end;quit" 


clear all;close all; parpool(24);
set(0,'DefaultAxesFontSize',14);

all_g_f = {
'cord.singlecontrolm0m1b.Anabaena_cylindrica_PCC_7122_uid183339.NC_019771.CCGCGG.data.txt',
'cord.singlecontrolm0m1b.Anabaena_cylindrica_PCC_7122_uid183339.NC_019771.CCTAGG.data.txt',
'cord.singlecontrolm0m1b.Anabaena_cylindrica_PCC_7122_uid183339.NC_019771.TACGTA.data.txt',
'cord.singlecontrolm0m1b.Anabaena_cylindrica_PCC_7122_uid183339.NC_019771.TCCGGA.data.txt',
'cord.singlecontrolm0m1b.Anabaena_cylindrica_PCC_7122_uid183339.NC_019771.TTCGAA.data.txt',
'cord.singlecontrolm0m1b.Chroococcidiopsis_thermalis_PCC_7203_uid183002.NC_019695.AATATT.data.txt',
'cord.singlecontrolm0m1b.Chroococcidiopsis_thermalis_PCC_7203_uid183002.NC_019695.CAATTG.data.txt',
'cord.singlecontrolm0m1b.Chroococcidiopsis_thermalis_PCC_7203_uid183002.NC_019695.GAATTC.data.txt',
'cord.singlecontrolm0m1b.Chroococcidiopsis_thermalis_PCC_7203_uid183002.NC_019695.TTGCAA.data.txt',
'cord.singlecontrolm0m1b.Chroococcidiopsis_thermalis_PCC_7203_uid183002.NC_019695.TTTAAA.data.txt',
'cord.singlecontrolm0m1b._Nostoc_azollae__0708_uid49725.NC_014248.AGGCCT.data.txt',
'cord.singlecontrolm0m1b._Nostoc_azollae__0708_uid49725.NC_014248.CCGCGG.data.txt',
'cord.singlecontrolm0m1b._Nostoc_azollae__0708_uid49725.NC_014248.GACGTC.data.txt',
'cord.singlecontrolm0m1b._Nostoc_azollae__0708_uid49725.NC_014248.TACGTA.data.txt',
'cord.singlecontrolm0m1b._Nostoc_azollae__0708_uid49725.NC_014248.TCCGGA.data.txt'
}



%parpool(4);

binSizeStart = 25;
binSizeEnd =7000;
binSizeInc = 25;
binCount = 200;
binSizes = binSizeStart:binSizeInc:binSizeEnd;

My_Rand_Seed = 33;
My_Rand_Seed_Start = 33;
num_rand = 200;

% all_GoF = zeros(length(all_g_f),length(binSizes))';
% all_Rp = zeros(length(all_g_f),length(binSizes))';
% all_GoF_c2 = zeros(length(all_g_f),length(binSizes))';
% all_Rp_c2 = zeros(length(all_g_f),length(binSizes))';
all_data = zeros(length(all_g_f),length(binSizes),9);
all_data_c2 = zeros(length(all_g_f),length(binSizes),9);
all_data_rand = zeros(length(all_g_f),length(binSizes),num_rand,9);

tic
%for j =[26:29,31:50]%,:length(all_g_f)
%for j =1:length(all_g_f)
%for j =11:length(all_g_f)
%for j =31:length(all_g_f)
%for j = [1:29,31:length(all_g_f)]
%for j = [39:length(all_g_f)]
%for j = [26:29,31:length(all_g_f)]
% ab 38, 39, 41
%for j = [42:length(all_g_f)]
%for j =41:length(all_g_f)
for j = 1:length(all_g_f)
    j
	file_cord = all_g_f{j};
	file_cord_c2 = strrep(file_cord,'m0','newC2');
	
	temp_loaded_file_cord = load(file_cord);
	g_size = temp_loaded_file_cord(1);
% 	if (g_size > 7000*binCount*2)
% 		binSizeEnd = min(7000,  floor(g_size/2/binCount));
% 		binSizes = binSizeStart:binSizeInc:binSizeEnd;
% 	else
% 		binSizeEnd = 7000;
% 		binSizes = binSizeStart:binSizeInc:binSizeEnd;
% 	end
	
	GoF_Array = zeros(length(binSizes),1);
	RPeriod_Array = zeros(length(binSizes),1);	
	data_815 = load(file_cord);
    N_815 = data_815(1);
    c_815 = data_815(2:end);
	c_GRand_allDist = getCircularDist(c_815,N_815);
	nM = size(c_815,1);
	
	%tic
	parfor i = 1:length(binSizes)
		if (i <= 200)
			continue;
		end
		%SameDMoutput = periodicity_matlab_v5sm_nSM_func_febGoF(file_cord,binCount,binSizes(i),8,0,0);
		SameDMoutput = periodicity_matlab_v5bsm_GRand_nSM_aprFunc(c_GRand_allDist,N_815,nM,binCount,binSizes(i));
		%GoF_Array(i) = SameDMoutput(7);
		%RPeriod_Array(i) = SameDMoutput(9);
        all_data(j,i,:) = SameDMoutput;
	end
	%toc
    %num_HIP1 = SameDMoutput(8);
    num_HIP1 = all_data(j,1,8);

	GoF_Array_rand = zeros(length(binSizes),1);
	RPeriod_Array_rand = zeros(length(binSizes),1);
	c_GRand_allDist = [];
	
	for jjj = 1:num_rand
		jjj
		rng(My_Rand_Seed_Start+jjj);
		c_GRand = sort(randsample(N_815,size(c_815,1)));
		c_GRand_allDist = getCircularDist(c_GRand,N_815);
		%tic
		parfor i = 1:length(binSizes)
			%SameDMoutput = periodicity_matlab_v5bsm_GRand_nSM_func_febGoF(c_GRand,N_815,binCount,binSizes(i),8,0,1);
			SameDMoutput = periodicity_matlab_v5bsm_GRand_nSM_aprFunc(c_GRand_allDist,N_815,nM,binCount,binSizes(i));
			all_data_rand(j,i,jjj,:) = SameDMoutput;
		end
		%toc
		%num_rand = num_HIP1;    
		%num_rand = all_data_rand(j,1,8);
		clear c_GRand_allDist;
		data_filename_temp = strcat('temp_matlab_save_6mer_july1.mat');
		save(data_filename_temp,'all_data','all_data_rand');
	end
	
  
	
	%fff
    toc
    %fff
    %data_filename = strcat('all_data.' ,'(',num2str(binCount),').mat');
	data_filename = strcat('all_data_6mer_july29.' ,'(',num2str(j),').mat');
    %save(data_filename);
	save(data_filename,'all_data','all_data_rand');
end



%ffff

%%% NO NEED TO DO BIN COUNT 100



