function Joptim=Balance_J(we,C)

Nnew=size(C,1);
dt=0.1;
tmax1=10000;
tspan1=0:dt:tmax1;

taon=100;
taog=10;
gamma=0.641;
sigma=0.01;
JN=0.15;
J=ones(Nnew,1);
I0=0.382;
Jexte=1.;
Jexti=0.7;
w=1.4;

curr=zeros(tmax1,Nnew);
delta=0.02*ones(Nnew,1);
for k=1:50000
 sn=0.001*ones(Nnew,1);
 sg=0.001*ones(Nnew,1);
 nn=1;
 j=0;
 for i=2:1:length(tspan1)
  xn=I0*Jexte+w*JN*sn+we*JN*C*sn-J.*sg;
  xg=I0*Jexti+JN*sn-sg;
  rn=phie(xn);
  rg=phii(xg);
  sn=sn+dt*(-sn/taon+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
  sn(sn>1) = 1;  
  sn(sn<0) = 0;      
  sg=sg+dt*(-sg/taog+rg./1000.)+sqrt(dt)*sigma*randn(Nnew,1);
  sg(sg>1) = 1;        
  sg(sg<0) = 0;
  j=j+1;
  if j==10
   curr(nn,:)=xn'-125/310;
   nn=nn+1;
   j=0;
  end
 end

 currm=mean(curr(1000:end,:),1);

 flag=0;
 for n=1:1:Nnew
  if abs(currm(n)+0.026)>0.005
   if currm(n)<-0.026 
    J(n)=J(n)-delta(n);
    delta(n)=delta(n)-0.001;
    if delta(n)<0.001;
       delta(n)=0.001;
    end
   else 
    J(n)=J(n)+delta(n);
   end
  else
   flag=flag+1;
  end
 end

 if flag==Nnew
  break;
 end
end

Joptim=J;
end
