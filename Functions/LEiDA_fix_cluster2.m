function [PTRANSITION,Pstates,LTime] = LEiDA_fix_cluster(BOLDsig,NumClusters,Center,TR)

[N_areas, Tmax]=size(BOLDsig);

% Preallocate variables to save FC patterns and associated information
Leading_Eig=zeros(Tmax,1*N_areas); % All leading eigenvectors

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

Phase_BOLD=zeros(N_areas,Tmax);

% Get the BOLD phase using the Hilbert transform
for seed=1:N_areas
    BOLDsig(seed,:)=BOLDsig(seed,:)-mean(BOLDsig(seed,:));
    signal_filt =filtfilt(bfilt,afilt,BOLDsig(seed,:));
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
    
    % Save V1 from all frames in all fMRI sessions in Leading eig
    t_all=t_all+1; % Update time
    Leading_Eig(t_all,:)=V1; %vertcat(V1,V2);
end

clear signal_filt iFC V1 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors
IDX=zeros(t_all,1);

for t=1:t_all
    for j=1:NumClusters
       di(j)=sqrt(sum((Leading_Eig(t,:)-Center(j,:)).^2));
%       di(j)=1-dot(Leading_Eig(t,:),Center(j,:))/norm(Leading_Eig(t,:))/norm(Center(j,:));
%        cc=1-corrcoef(Leading_Eig(t,:),Center(j,:));
%        di(j)=cc(2);
    end
    [aux indmin]=min(di);
    IDX(t)=indmin;
end

Pstates=zeros(1,NumClusters);
for c=1:NumClusters
    Pstates(c)=mean(IDX==c);
    % Mean Lifetime
    Ctime_bin=IDX==c;
    
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
    LTime(c)=mean(C_Durations)*TR;
end

Pstates=Pstates/sum(Pstates);

[~, ind_sort]=sort(Pstates,'descend');

PTRANSITION=zeros(NumClusters,NumClusters);
i=1;
for c1=ind_sort
    j=1;
    for c2=ind_sort
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1)));
        end
        j=j+1;
    end
    i=i+1;
end

 