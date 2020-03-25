clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%


Wexc=0.:0.025:0.5;
Winh=0.:0.025:0.5;
for s=1:441
    fileName = sprintf('Gw347nr_%03d.mat',s);
    load(fileName);
    [IWexc IWinh]=ind2sub([length(Wexc),length(Winh)],s);
    fittpla_all(IWexc,IWinh)=fittpla;
    klpstatespla_all(IWexc,IWinh)=klpstatespla;
    entropydistpla_all(IWexc,IWinh)=entropydistpla;
    errorlifetimepla_all(IWexc,IWinh)=errorlifetimepla;
    fittlsd_all(IWexc,IWinh)=fittlsd;
    klpstateslsd_all(IWexc,IWinh)=klpstateslsd;
    entropydistlsd_all(IWexc,IWinh)=entropydistlsd;
    errorlifetimelsd_all(IWexc,IWinh)=errorlifetimelsd;
end

save neuromodulator_2D_Psilo347nr.mat Wexc Winh fittpla_all klpstatespla_all entropydistpla_all errorlifetimepla_all fittlsd_all klpstateslsd_all entropydistlsd_all errorlifetimelsd_all;


figure
imagesc(Winh,Wexc,klpstateslsd_all)

figure
imagesc(Winh,Wexc,errorlifetimelsd_all)

figure
subplot(1,2,1);
imagesc(Winh,Wexc,klpstateslsd_all)
subplot(1,2,2);
plot(sort(klpstateslsd_all(:),'descend'))
ylim([0 0.031])

figure
subplot(1,2,1);
imagesc(Winh,Wexc,errorlifetimelsd_all)
subplot(1,2,2);
plot(sort(errorlifetimelsd_all(:),'descend'))
ylim([0 7.74])