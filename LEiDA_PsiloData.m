clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%

addpath('/Users/josephine/Documents/MATLAB/Neurotrans_2/Functions')

load psilotc.mat;

NSUB=9; % Number of Subjects

%CASE 3= placebo ; CASE 4= Psilo

for nsub=1:NSUB
    CASE=3;
    tc_aal_new{nsub,1}=LR_version_symm(tc_aal{nsub,CASE});
    CASE=4;
    tc_aal_new{nsub,2}=LR_version_symm(tc_aal{nsub,CASE});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 - Compute the Leading Eigenvectors from the BOLD datasets
disp('Processing the eigenvectors from BOLD data')
% Load here the BOLD data (which may be in different formats)
% Here the BOLD time courses in AAL parcellation are organized as cells,
% where tc_aal{1,1} corresponds to the BOLD data from subject 1 in
% condition 1 and contains a matrix with lines=N_areas and columns=Tmax.

[n_Subjects, n_cond]=size(tc_aal_new);
[N_areas, Tmax]=size(tc_aal_new{1,1});

% Parameters of the data
TR=3.;  % Repetition Time (seconds)

% Preallocate variables to save FC patterns and associated information
Leading_Eig=zeros(Tmax*n_Subjects,1*N_areas); % All leading eigenvectors
Time_all=zeros(2, n_Subjects*Tmax); % vector with subject nr and cond at each t
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

for s=1:n_Subjects
    for cond=1:n_cond
        [N_areas, Tmax]=size(tc_aal_new{s,cond});
        % Get the BOLD signals from this subject in this condition
        BOLD = tc_aal_new{s,cond};
        Phase_BOLD=zeros(N_areas,Tmax);
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            ts=demean(detrend(BOLD(seed,:)));
            signal_filt =filtfilt(bfilt,afilt,ts);
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        for t=1:Tmax
            
            %Calculate the Instantaneous FC (BOLD Phase Synchrony)
            iFC=zeros(N_areas);
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            % Get the leading eigenvector
            
            [V1,~]=eigs(iFC,1);
            % Make sure the largest component is negative
            if mean(V1>0)>.5
                V1=-V1;
            elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
                V1=-V1;
            end           
            % Save V1 from all frames in all fMRI sessions in Leading eig
            t_all=t_all+1; % Update time
            Leading_Eig(t_all,:)=V1; %vertcat(V1,V2);
            Time_all(:,t_all)=[s cond]; % Information that at t_all, V1 corresponds to subject s in a given condition
        end
    end
end
clear BOLD tc_aal signal_filt iFC VV V1 V2 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors

disp('Clustering the eigenvectors into')
% Leading_Eig is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: Tmax*n_Subjects are all time points (independent observations)

% Set maximum/minimum number of clusters
% There is no fixed number of states the brain can display
% Extending depending on the hypothesis of each work
maxk=3;
mink=3;
rangeK=mink:maxk;

% Set the parameters for Kmeans clustering
Kmeans_results=cell(size(rangeK));
for k=1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' clusters'])
    [IDX, C, SUMD, D]=kmeans(Leading_Eig,rangeK(k),'Distance','sqeuclidean','Replicates',200,'Display','off');
    Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric collumn vectos
    Kmeans_results{k}.C=C;       % Cluster centroids (FC patterns)
    Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D;       % Distance from each point to every centroid
    ss=silhouette(Leading_Eig,IDX,'sqeuclidean');
    sil(k)=mean(ss);
end

distM_fcd=squareform(pdist(Leading_Eig,'euclidean'));
dunn_score=zeros(length(rangeK),1);
for j=1:length(rangeK)
    dunn_score(j)=dunns(rangeK(j),distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(rangeK(j)) ' clusters'])
end

[~,ind_max]=max(dunn_score);

disp(['Best clustering solution: ' num2str(rangeK(ind_max)) ' clusters']);

%% 3 - Analyse the Clustering results

% For every fMRI scan calculate probability and lifetimes of each pattern c.
P=zeros(n_cond,n_Subjects,maxk-mink+1,maxk);
LT=zeros(n_cond,n_Subjects,maxk-mink+1,maxk);

for k=1:length(rangeK)
    for cond=1:n_cond   
        for s=1:n_Subjects
            
            % Select the time points representing this subject and task
            T=((Time_all(1,:)==s)+(Time_all(2,:)==cond))>1;
            Ctime=Kmeans_results{k}.IDX(T);
            
            for c=1:rangeK(k)
                % Probability
                P(cond,s,k,c)=mean(Ctime==c);
                
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(cond,s,k,c)=mean(C_Durations)*TR;
            end
        end
    end
end

P_pval=zeros(maxk-mink+1,maxk);
LT_pval=zeros(maxk-mink+1,maxk);

disp('Test significance between Placebo and LSD')
for k=1:length(rangeK)
    disp(['Now running for ' num2str(k) ' clusters'])
    for c=1:rangeK(k)
        % Compare Probabilities
        a=squeeze(P(1,:,k,c));  % Vector containing Prob of c in Baselineline
        b=squeeze(P(2,:,k,c));  % Vector containing Prob of c in LSD
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        P_pval(k,c)=min(stats.pvals);
        
        % Comapre Lifetimes
        a=squeeze(LT(1,:,k,c));  % Vector containing Prob of c in Baselineline
        b=squeeze(LT(2,:,k,c));  % Vector containing Prob of c in LSD
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        LT_pval(k,c)=min(stats.pvals);
    end
end
disp('%%%%% LEiDA SUCCESSFULLY COMPLETED %%%%%%%')

%% 4 - Plot FC patterns and stastistics between groups

disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])

for k=1:length(rangeK)
    Pmin_pval(k)=min(P_pval(k,1:rangeK(k)));
end

[pmin_pval k]=min(Pmin_pval);
disp(['Note: The minimal K with significant difference is detected with K=' num2str(rangeK(k)) ' (p=' num2str(pmin_pval) ')'])

numsig=-1;
for k=1:length(rangeK)
    nums=length(find(P_pval(k,1:rangeK(k))<0.05));
    if nums>numsig
        numsig=nums;
    end
end
disp(['Note: The K with max number of significant difference is detected with K=' num2str(rangeK(k)) ' (n=' num2str(numsig) ')'])

disp(['Note: Silohuette optimum:'])
sil

%%% Graphics

K = input('Number of clusters: ');
Number_Clusters=K;
Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Clusters are sorted according to their probability of occurrence
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend');

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);
% Order=[1:2:N N:-2:2];

Vemp=V;
P1emp=squeeze(P(1,:,k,ind_sort));
P2emp=squeeze(P(2,:,k,ind_sort));
LT1emp=squeeze(LT(1,:,k,ind_sort));
LT2emp=squeeze(LT(2,:,k,ind_sort));

meanp1emp=squeeze(mean(P1emp));
meanp2emp=squeeze(mean(P2emp));
[~, ind_sort1]=sort(meanp1emp,'descend');
[~, ind_sort2]=sort(meanp2emp,'descend');


PTR1emp=zeros(K);
PTR2emp=zeros(K);
n_sub1=zeros(K,1);
n_sub2=zeros(K,1);
for psylo=1:2 
    for s=1:n_Subjects
        % Select the time points representing this subject and LSD
        T=((Time_all(1,:)==s)+(Time_all(2,:)==psylo))>1;
        Ctime=Kmeans_results{k}.IDX(T);
        
        if psylo==1
            i=1;
            for c1=ind_sort1
                j=1;
                for c2=ind_sort1
                    sumatr=0;
                    for t=1:length(Ctime)-1
                        if Ctime(t)==c1 && Ctime(t+1)==c2
                            sumatr=sumatr+1;
                        end
                    end
                    if length(find(Ctime(1:length(Ctime)-1)==c1)) ~= 0
                        PTR1emp(i,j)=PTR1emp(i,j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
                    end
                    j=j+1;
                end
                if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                    n_sub1(i)=n_sub1(i)+1;
                end
                i=i+1;
            end
        end
        
        if psylo==2
            i=1;
            for c1=ind_sort2
                j=1;
                for c2=ind_sort2
                    sumatr=0;
                    for t=1:length(Ctime)-1
                        if Ctime(t)==c1 && Ctime(t+1)==c2
                            sumatr=sumatr+1;
                        end
                    end
                    if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                        PTR2emp(i,j)=PTR2emp(i,j)+sumatr/(length(find(Ctime(1:length(Ctime)-1)==c1)));
                    end
                    j=j+1;
                end
                if length(find(Ctime(1:length(Ctime)-1)==c1)) ~=0
                    n_sub2(i)=n_sub2(i)+1;
                end
                
                i=i+1;
            end
        end
       
    end
end
for i=1:K
    PTR1emp(i,:)=PTR1emp(i,:)/n_sub1(i);
    PTR2emp(i,:)=PTR2emp(i,:)/n_sub2(i);
end
    
save empiricalLEiDA_Psilo3_0407.mat Vemp P1emp P2emp LT1emp LT2emp PTR1emp PTR2emp Number_Clusters;

