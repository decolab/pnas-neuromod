clc; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%

load  empiricalLEiDA_Psilo3_0407.mat;
load  raphesymnorm2.mat;


P1emp=mean(P1emp);
P2emp=mean(P2emp);

LT1emp=mean(LT1emp);
LT2emp=mean(LT2emp);

load all_SC_FC_TC_76_90_116.mat;
load mean5T_all;

C=sc90;
C=C/max(max(C))*0.2;
N=90;

Receptor=symm_mean5HT2A/max(symm_mean5HT2A);
rng(1);
index=randperm(90);
Receptor1=Receptor(index);

for s=1
    rng(s);
    
    TR=3;
    NumClusters=Number_Clusters;
    
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
    
    Jlsd=0.1;
    
    Tmax=2000;
    boldstep=TR*1000;
    
    %%%%%%%%%%%%
    %% Optimize
    %%
    tlsd=120;
    we=1.6;
    wexc=0.3;
    winh=0.1;
    display('Dynamics out');
    
    J=Balance_J(we,C);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    
    j=1;
    for t=0:dt:(1000*100*TR)
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
        Ylsd(:,j)=Ilsd;
        j=j+1;
    end
    
    Ilsd=mean(Ylsd(:,end-10000:end),2);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    %         Ilsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('dynout_%03d.mat',s),'klpstateslsd2','errorlifetimelsd2','entropydistlsd2');
    
    
    %%%%%%%%%%%%
    display('Dynamics in');
    
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2dyn=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2dyn=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2dyn=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('dynin_%03d.mat',s),'klpstateslsd2dyn','errorlifetimelsd2dyn','entropydistlsd2dyn');
    
    
    %%%%%%%%%%%%
    
    display('Random');
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor1.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor1.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2rnd=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2rnd=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2rnd=EntropyMarkov(PTR2emp,PTRsim);
    P2emp
    Pstates
    Ind
    
    
    save(sprintf('reshuf_%03d.mat',s),'klpstateslsd2rnd','errorlifetimelsd2rnd','entropydistlsd2rnd');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('Receptor 1A');
    
    Receptor1=symm_mean5HT1A/max(symm_mean5HT2A);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor1.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor1.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd21A=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd21A=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd21A=EntropyMarkov(PTR2emp,PTRsim);
    
    
    
    save(sprintf('rec1A_%03d.mat',s),'klpstateslsd21A','errorlifetimelsd21A','entropydistlsd21A');
    
    %%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('Receptor 1B');
    
    Receptor1=symm_mean5HT1B/max(symm_mean5HT2A);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor1.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor1.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd21B=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd21B=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd21B=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('rec1B_%03d.mat',s),'klpstateslsd21B','errorlifetimelsd21B','entropydistlsd21B');
    
    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('Receptor T4');
    
    Receptor1=symm_mean5HT4/max(symm_mean5HT2A);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor1.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor1.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2T4=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2T4=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2T4=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('recT4_%03d.mat',s),'klpstateslsd2T4','errorlifetimelsd2T4','entropydistlsd2T4');
    
    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display('Receptor TT');
    
    Receptor1=symm_mean5HTT/max(symm_mean5HT2A);
    
    neuro_act=zeros(round(1000*(Tmax-1)*TR+1),N);
    sn=0.001*ones(N,1);
    sg=0.001*ones(N,1);
    rn=phie(I0*Jexte+w*JN*sn-J.*sg);
    Ilsd=zeros(N,1);
    ylsd=zeros(N,1);
    nn=1;
    for t=0:dt:(1000*(Tmax-1)*TR)
        xn=I0*Jexte+w*JN*sn+we*JN*C*sn+wexc*Receptor1.*Ilsd-J.*sg;
        xg=I0*Jexti+JN*sn+winh*Receptor1.*Ilsd-sg;
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2TT=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2TT=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2TT=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('recTT_%03d.mat',s),'klpstateslsd2TT','errorlifetimelsd2TT','entropydistlsd2TT');
    
    %%%%%%%%%%%%
    display('Dynamics in we wi 0');
    wexc=0;
    winh=0;
    
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
    
    bds=BOLD_act(5*boldstep:boldstep:end,:);
    %%%
    
    
    %%%% KL dist between PTR2emp
    
    [PTRsim,Pstates,LTime]=LEiDA_fix_cluster2(bds',NumClusters,Vemp,TR);
    
    errorlifetimelsd2dyn0=sqrt(sum((LT2emp-LTime).^2)/length(LTime));
    Ind=find(Pstates~=0);
    klpstateslsd2dyn0=0.5*(sum(Pstates(Ind).*log(Pstates(Ind)./P2emp(Ind)))+sum(P2emp(Ind).*log(P2emp(Ind)./Pstates(Ind))));
    entropydistlsd2dyn0=EntropyMarkov(PTR2emp,PTRsim);
    
    
    save(sprintf('zero_%03d.mat',s),'klpstateslsd2dyn0','errorlifetimelsd2dyn0','entropydistlsd2dyn0');
    
end
