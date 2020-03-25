clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%


for s=1:9
    fileName = sprintf('w00_%03d.mat',s);
    load(fileName);
    Pstates_all00(s,:)=Pstates;
    fileName = sprintf('wopt_%03d.mat',s);
    load(fileName);
    Pstates_allopt(s,:)=Pstates;   
end

figure
bar(mean(Pstates_all00))
hold
errorbar(mean(Pstates_all00),std(Pstates_all00)/sqrt(9),'LineStyle','none')

figure
bar(mean(Pstates_allopt))
hold
errorbar(mean(Pstates_allopt),std(Pstates_allopt)/sqrt(9),'LineStyle','none')