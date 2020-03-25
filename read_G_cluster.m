clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Dynamic coupling of Whole-Brain Neuronal and Neurotransmitter Systems
%     Kringelbach, M. L., Cruzat, J., Cabral, J., Knudsen, G. M.,
%       Carhart-Harris, R. L., Whybrow, P. C., Logothetis N. K. & Deco, G.
%         (2020) Proceedings of the National Academy of Sciences

%   Barcelona?Spain, March, 2020.

%%%%%%

WE=0.5:0.01:2;

for s=1:151
    fileName = sprintf('Gp_%03d.mat',s);
    load(fileName);
    fitt_all(s)=fitt;
    klpstatesplacebo_all(s)=klpstatesplacebo;
    entropydistplacebo_all(s)=entropydistplacebo;
    errorlifetimeplacebo_all(s)=errorlifetimeplacebo;
    Coptim(s,:,:)=Cnew;
end

save optimizedmfplacebo_Psilo5.mat WE fitt_all klpstatesplacebo_all entropydistplacebo_all errorlifetimeplacebo_all;

figure
subplot(3,1,1);
plot(WE,fitt_all,'b');
subplot(3,1,2);
plot(WE,klpstatesplacebo_all,'k');
subplot(3,1,3);
plot(WE,errorlifetimeplacebo_all,'r');



