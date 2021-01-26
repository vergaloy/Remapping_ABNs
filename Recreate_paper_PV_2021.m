%% Prepare Data
%hyp: hypnogram 0=W, 0.5=NREM, 1=REM.
bin=1;  %bin size
sf=5;   %frame rate
[mice_sleep,N]=separate_data_by_conditions(data); %HC,preS,postS,Consolidation,Test
[D,a,at]=bin_mice_sleep(mice_sleep,bin); % Bin data (bin size=1s)
%% Figure 1
[increasing,decreasing]=ABNs_responding_to_preS_and_postS(D,10);
%% Figure 2
clus=cluster_activity_vectors(at); 
%% Figure 3
out=get_ABNs_remapping_in_time(a,D,N,hyp,sf,bin); 