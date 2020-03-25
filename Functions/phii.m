function result=phie(x)
g=0.087;
I=177.;
c=615.; 
y=c*x-I;
 if y~=0
  result = y./(1-exp(-g*y));
 else
  result=0;
 end
end
