clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%

load  empiricalLEiDA_Psilo3_0407.mat;

P1emp=mean(P1emp);
LT1emp=mean(LT1emp);
P2emp=mean(P2emp);
LT2emp=mean(LT2emp);

load all_SC_FC_TC_76_90_116.mat;
load mean5HT2A_bindingaal.mat

C=sc90;
C=C/max(max(C))*0.2;
N=90;

Receptor=mean5HT2A_aalsymm(:,1)/max(mean5HT2A_aalsymm(:,1));

load psilotc.mat;

NSUB=9;
Isubdiag = find(tril(ones(N),-1));

%%Here is my example of my data (are fdifferent cases...
CASE=3; %% 3 placebo ; 4 Psilo
TC1=LR_version_symm(tc_aal{1,CASE});
TC2=LR_version_symm(tc_aal{2,CASE});
TC3=LR_version_symm(tc_aal{3,CASE});
TC4=LR_version_symm(tc_aal{4,CASE});
TC5=LR_version_symm(tc_aal{5,CASE});
TC6=LR_version_symm(tc_aal{6,CASE});
TC7=LR_version_symm(tc_aal{7,CASE});
TC8=LR_version_symm(tc_aal{8,CASE});
TC9=LR_version_symm(tc_aal{9,CASE});

xs=eval(sprintf('TC%d',1));
Tmax=size(xs,2);

TR=3.;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = 0.04;           % lowpass frequency of filter
fhi = 0.07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter


%%%%%%%%%%%%%%

for nsub=1:NSUB
    signaldata=eval(sprintf('TC%d', nsub));
    FCemp2(nsub,:,:)=corrcoef(signaldata');
    Phase_BOLD_data=zeros(N,Tmax);
    for seed=1:N
        signaldata(seed,:)=signaldata(seed,:)-mean(signaldata(seed,:));
        signal_filt_data =filtfilt(bfilt2,afilt2,signaldata(seed,:));
        Phase_BOLD_data(seed,:) = angle(hilbert(signal_filt_data));
    end
    
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end
    FCphasesemp2(nsub,:,:)=squeeze(mean(iFC));
end
FCphasesemp=squeeze(mean(FCphasesemp2));
FCemp=squeeze(mean(FCemp2));

%%%%%%%%%%%%%%%%%%
dtt   = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
dt=0.1;

taon=100;
taog=10;
gamma=0.641;
sigma=0.01;
JN=0.15;
I0=0.382;
Jexte=1.;
Jexti=0.7;
w=1.4;

Tmax=2000;
boldstep=TR*1000;

%%%%%%%%%%%%
%% Optimize
%%
iwe=1;
WE=0.5:0.01:2;

for s = 1:(size(WE,2))
    we=WE(s);
    
    Cnew=C;
    
    J=Balance_J(we,Cnew);
    
    
    for iter=1:10
        neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
        sn=0.001*ones(N,1);
        sg=0.001*ones(N,1);
        nn=1;
        for t=0:dt:(1000*(Tmax-1)*TR)
            xn=I0*Jexte+w*JN*sn+we*JN*Cnew*sn-J.*sg;
            xg=I0*Jexti+JN*sn-sg;
            rn=phie(xn);
            rg=phii(xg);
            sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(N,1);
            sn(sn>1) = 1;
            sn(sn<0) = 0;
            sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(N,1);
            sg(sg>1) = 1;
            sg(sg<0) = 0;
            j=j+1;
            if abs(mod(t,1))<0.01
                neuro_act(nn,:)=rn';
                nn=nn+1;
            end
        end
        nn=nn-1;
        
        %%%% BOLD empirical
        % Friston BALLOON MODEL
        T = nn*dtt; % Total time in seconds
        
        B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
        BOLD_act = zeros(length(B),N);
        BOLD_act(:,1) = B;
        
        for nnew=2:N
            B = BOLD(T,neuro_act(1:nn,nnew));
            BOLD_act(:,nnew) = B;
        end
        
        bds=BOLD_act(boldstep:boldstep:end,:);
        %%%
        
        FC_simul=corrcoef(bds);
        cc=corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)));
        fitt2(iwe)=cc(2);
        
        %%%% KL dist between PTR2emp
        
        [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
        
        errorlifetimeplacebo2(iwe)=sqrt(sum((LT1emp-LTime).^2)/length(LTime))
        
        klpstatesplacebo2(iwe)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)))
        
        entropydistplacebo2(iwe)=EntropyMarkov2(PTR1emp,PTRsim,P1emp,Pstates)
        iwe=iwe+1;
        
    end
    
    fitt=mean(fitt2);
    errorlifetimeplacebo=mean(errorlifetimeplacebo2);
    klpstatesplacebo=mean(klpstatesplacebo2);
    entropydistplacebo=mean(entropydistplacebo2);
    
    
    save(sprintf('Gp_%03d.mat',s),'fitt','errorlifetimeplacebo','klpstatesplacebo','entropydistplacebo','Cnew');
    
end