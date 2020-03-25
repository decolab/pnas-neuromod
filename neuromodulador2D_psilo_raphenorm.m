clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%

load empiricalLEiDA_Psilo3_0407.mat;
load raphesymnorm2.mat;

P1emp=mean(P1emp);
P2emp=mean(P2emp);

LT1emp=mean(LT1emp);
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


%%%%%%%%%%%%%%

for nsub=1:NSUB
    signaldata=eval(sprintf('TC%d', nsub));
    FCemp2(nsub,:,:)=corrcoef(signaldata');
end
FCemppla=squeeze(mean(FCemp2));

%%%%%%%%%%%%%%%
CASE=4; %% 3 placebo ; 4 Psilo
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

for nsub=1:NSUB
    signaldata=eval(sprintf('TC%d', nsub));
    FCemp2(nsub,:,:)=corrcoef(signaldata');
end
FCemplsd=squeeze(mean(FCemp2));

%%%%%%%%%%%%%%%%%%
dtt = 1e-3;   % Sampling rate of simulated neuronal activity (seconds)
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

Jlsd=0.1;

Tmax=2000;
boldstep=TR*1000;

%%%%%%%%%%%%
%% Optimize
%%
Wexc=0.:0.025:0.5;
Winh=0.:0.025:0.5;

tlsd=120;


for s = 1:(length(Wexc)*length(Winh))
    [IWexc IWinh]=ind2sub([length(Wexc),length(Winh)],s);
    we=1.6;             % This is the minimun value found when optimized for the placebo
    wexc=Wexc(IWexc);
    winh=Winh(IWinh);
    
    J=Balance_J(we,C);
    for iter=1:10
        neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
        sn=0.001*ones(N,1);
        sg=0.001*ones(N,1);
        rn=phie(I0*Jexte+w*JN*sn-J.*sg);
        Ilsd=zeros(N,1);
        ylsd=zeros(N,1);
        nn=1;
        for t=0:dt:(1000*(Tmax-1)*TR)
            xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor.*Ilsd-J.*sg;
            xg=I0*Jexti+JN*sn+winh*Receptor.*Ilsd-sg;
            ylsd=ylsd+dt*(5*raphe*raphe'*rn-1300*ylsd./(170+ylsd));
            Ilsd=Ilsd+dt/tlsd*(-Ilsd+Jlsd./(1+exp((-log10(ylsd)+1)/0.1)));
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
        cc=corrcoef(atanh(FCemppla(Isubdiag)),atanh(FC_simul(Isubdiag)));
        fittpla2(iter)=cc(2);
        cc=corrcoef(atanh(FCemplsd(Isubdiag)),atanh(FC_simul(Isubdiag)));
        fittlsd2(iter)=cc(2);
        
        %%%% KL dist between PTR2emp
        
        [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
        
        errorlifetimepla2(iter)=sqrt(sum((LT1emp-LTime).^2)/length(LTime));
        %% klpstatespla2(iter)=1/sqrt(2)*sqrt(sum((sqrt(Pstates)-sqrt(P1emp)).^2));
        klpstatespla2(iter)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)))
        
        entropydistpla2(iter)=EntropyMarkov(PTR1emp,PTRsim);
        
        errorlifetimelsd2(iter)=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
        %% klpstateslsd2(iter)=1/sqrt(2)*sqrt(sum((sqrt(Pstates)-sqrt(P2emp)).^2));
        klpstateslsd2(iter)=0.5*(sum(Pstates.*log(Pstates./P2emp))+sum(P2emp.*log(P2emp./Pstates)))
        
        entropydistlsd2(iter)=EntropyMarkov(PTR2emp,PTRsim);
        
    end
    
    fittpla=mean(fittpla2);
    fittlsd=mean(fittlsd2);
    errorlifetimepla=mean(errorlifetimepla2);
    klpstatespla=mean(klpstatespla2);
    entropydistpla=mean(entropydistpla2);
    errorlifetimelsd=mean(errorlifetimelsd2);
    klpstateslsd=mean(klpstateslsd2);
    entropydistlsd=mean(entropydistlsd2);
    
    
    save(sprintf('Gw347nr_%03d.mat',s),'fittpla','fittlsd','errorlifetimepla','klpstatespla','entropydistpla','errorlifetimelsd','klpstateslsd','entropydistlsd');
    
end