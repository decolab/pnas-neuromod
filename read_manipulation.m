clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%


for s=1:100
    fileName = sprintf('dynout_%03d.mat',s);
    load(fileName);
    klpstateslsd2_all(s)=klpstateslsd2;
    errorlifetimelsd2_all(s)=errorlifetimelsd2;
    fileName = sprintf('dynin_%03d.mat',s);
    load(fileName);
    klpstateslsd2dyn_all(s)=klpstateslsd2dyn;
    errorlifetimelsd2dyn_all(s)=errorlifetimelsd2dyn;
    fileName = sprintf('reshuf_%03d.mat',s);
    load(fileName);
    klpstateslsd2rnd_all(s)=klpstateslsd2rnd;
    errorlifetimelsd2rnd_all(s)=errorlifetimelsd2rnd;
    fileName = sprintf('zero_%03d.mat',s);
    load(fileName);
    klpstateslsd2dyn0_all(s)=klpstateslsd2dyn0;
    errorlifetimelsd2dyn0_all(s)=errorlifetimelsd2dyn0;
    fileName = sprintf('rec1A_%03d.mat',s);
    load(fileName);
    klpstateslsd21A_all(s)=klpstateslsd21A;
    errorlifetimelsd21A_all(s)=errorlifetimelsd21A;
    fileName = sprintf('rec1B_%03d.mat',s);
    load(fileName);
    klpstateslsd21B_all(s)=klpstateslsd21B;
    errorlifetimelsd21B_all(s)=errorlifetimelsd21B;
    fileName = sprintf('recTT_%03d.mat',s);
    load(fileName);
    klpstateslsd2TT_all(s)=klpstateslsd2TT;
    errorlifetimelsd2TT_all(s)=errorlifetimelsd2TT;
    fileName = sprintf('recT4_%03d.mat',s);
    load(fileName);
    klpstateslsd2T4_all(s)=klpstateslsd2T4;
    errorlifetimelsd2T4_all(s)=errorlifetimelsd2T4;
end

figure
boxplot([klpstateslsd2dyn_all' klpstateslsd2dyn0_all' klpstateslsd2_all' klpstateslsd2rnd_all' klpstateslsd21A_all' klpstateslsd21B_all' klpstateslsd2T4_all' klpstateslsd2TT_all']);
figure
boxplot([errorlifetimelsd2dyn_all' errorlifetimelsd2dyn0_all' errorlifetimelsd2_all' errorlifetimelsd2rnd_all' errorlifetimelsd21A_all' errorlifetimelsd21B_all' errorlifetimelsd2T4_all' errorlifetimelsd2TT_all']);

        a=klpstateslsd2dyn0_all; 
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pkl0=min(stats.pvals);
        a=errorlifetimelsd2dyn0_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptime0=min(stats.pvals);
    
        a=klpstateslsd2_all; 
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pkl=min(stats.pvals);
        a=errorlifetimelsd2_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptime=min(stats.pvals);
        
        a=klpstateslsd2rnd_all; 
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pklrnd=min(stats.pvals);
        a=errorlifetimelsd2rnd_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptimernd=min(stats.pvals);

        a=klpstateslsd21A_all; 
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pklr1A=min(stats.pvals);
        a=errorlifetimelsd21A_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptime1A=min(stats.pvals);
        
        a=klpstateslsd21B_all;
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pklr1B=min(stats.pvals);
        a=errorlifetimelsd21B_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptime1B=min(stats.pvals);

        a=klpstateslsd2T4_all;
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pklrT4=min(stats.pvals);
        a=errorlifetimelsd2T4_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptimeT4=min(stats.pvals);
        
        a=klpstateslsd2TT_all;
        b=klpstateslsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        pklrTT=min(stats.pvals);
        a=errorlifetimelsd2TT_all; 
        b=errorlifetimelsd2dyn_all; 
        stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
        ptimeTT=min(stats.pvals);
        pkl0
        ptime0
        pkl
        ptime
        pklrnd
        ptimernd
        pklr1A
        ptime1A
        pklr1B
        ptime1B
        pklrT4
        ptimeT4
        pklrTT
        ptimeTT
        
save manipulation.mat;
