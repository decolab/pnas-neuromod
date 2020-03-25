function result=phie(x)
g=0.16;
I=125.;
c=310.;
y=c*x-I;
 if y~=0
  result = y./(1-exp(-g*y));
 else
  result=0;
 end
end
