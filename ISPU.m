%................explaination for Input parameters of ISPU......
%y is  vector of the desire signal; x is vector of the input signal;
%h is IR of the unknown system. It is not necessary for the ISPU. It is
%only used for caculating the misalignment.
%mu is the step size, M is the number of blocks.
%p is the update probability of the inactive blocks
%wk is the initialed vector of the adaptive filter;

%.....explaination for output of ISPU...........
%mis is vector of misalignment
%filter is the final state of the adaptive filter
%r1 is the vector keeps the curves of the ratio between the number of addition required by SPU and NLMS 
%r2 is the vector keeps the curves of the ratio between the number of comparison required by SPU and the length of filter.

%...................An example for running the ISPU...............
% h=[zeros(1,128) randn(1,64) zeros(1,64)];
% xk=10*randn(1,10000);
% xk=filter(1,[1 -0.7],xk);%Gererate a colored noise;
% y=conv(xk,h);
% noise=(sqrt(var(y)/1000))*randn(1,length(y));%SNR=30dB
% echo=noise+y;
% [mis,r1,r2,filter]=ISPU(echo,xk,h,0.5,32,1/8,zeros(1,256));


function  [mis,r1,r2,filter]=ISPU(y,x,h,mu,M,p,wk)
%................initialization for SPU........................
h=MyReshape(h);
wk=MyReshape(wk);
global w;
w=wk;

N=length(h);
L=N/M;
len=min(length(x),length(y));
xk=x(N:-1:1);
v1=5*var(x);
v1=v1+sum(xk.*xk); %denominator for the adaptation equation of NLMS ;
global l2n;%vector to keep the l2 norm of each input blocks;
l2n=L2nofEachBlock(M,L,xk);

yk=0;
ek=0;
den=sum(h.*h); %denominator for the Misalignment.
mis=zeros(1,len-N-1); % vector  to keep the misalignment of algorithm in each adaptation;
jjj=1;%just and index for keeping the misalignment in the vector mis;
global Ai;  %Vector to keep the indexes of the located active blocks;
Ai=zeros(1,M);
Alen=0;  %the number of the located active blocks;
global Zi;  %Vector to keep the indexes of the located zero blocks;
Zi=zeros(1,M);
Zlen=0;   %the number of the located active blocks;
nz=0;  %nz is the number of the selected zero blocks in each adaptation;
nu=0;  %nu is the number of the selected adapting blocks, i.e., including all the active blocks and the selected zero blocks,  in each adaptation;
out_x=xk(N);
thre=M/2;
%thre=floor(M*(1-1/L)-1.5);%

flag_NLMS=0;
flag_SPU=0;
if all (wk==0)
    flag_NLMS=100;
    flag_SPU=0;     
else
    flag_SPU=100;
    flag_NLMS=0;
    l2nwk=L2nofEachBlock(M,L,wk);
    [minvalue,minindex]=min(l2nwk);
    Zlen=1;
    Zi(Zlen)=minindex; 
    [Alen,Zlen]=LocateActiveBlock(Zlen,L,M);    
    nz=Calculate_nz(p,Zlen);  
end
j=1;
r1=mis;
r2=r1;
%................initialization for SPU........................

%[Alen,Zlen,nz]=InitialActiveBlock(M,p);
%j=1;
for index=N+1:len  
     xk=x(index:-1:index-N+1); %xk is the input vector of the adaptive algorithm.
     v1=v1-out_x*out_x+xk(1)*xk(1); %v1 is the denominator for NLMS. 
     maintainL2n(M,L,out_x,xk);%maintian the L2 norm of each input block by subtracting off the outgoing sample and adding the incoming sample
     if (flag_SPU>0)         
         nu=GenerateAdaptingBlocksIndex(j,Alen,Zlen,nz); 
         yk=GenerateOutput(nu,xk,L);% generate the filter output of SPU;
         v2=GenerateDenominatorSPU(nu);   %generate the denominator for the SPU adaptation ;
         ek=y(index)-yk;
         UpdateCoefficients(mu,ek,v2,xk,nu,L);   
         mis(jjj)=CalculateMisalignment(nu,h,L,den);%Misalignment is calculated to evalute the algorithm performance. It is not necessary to the algorithm application.
         r1(jjj)=Alen+2*nz;
         r2(jjj)=0;
         j=j+nz; 
         if((j>Zlen))
             [Alen,Zlen,n]=LocateActiveBlock(Zlen,L,M); 
             nz=Calculate_nz(p,Zlen);
             j=1;
             r2(jjj)=n;
         end 
         flag_SPU=flag_SPU-1;
     end
     if (flag_NLMS>0)
        yk=sum(xk.*w);
        ek=y(index)-yk;
        temp1=mu*ek/v1;
        w=w+temp1.*xk;
        temp2=w-h;
        m=sum(temp2.*temp2);
        mis(jjj)=10*log10(m/den);
        r1(jjj)=M;
        r2(jjj)=0;
        flag_NLMS=flag_NLMS-1;        
     end      
     if ((flag_SPU==0)&&(flag_NLMS==0))
         [Alen,Zlen,nz,flag_NLMS]=InitialActiveBlock(M,L,thre,Zlen,p,Alen,nz);
         if (flag_NLMS==100)
             flag_SPU=0;
         else
             flag_SPU=100;
         end         
     end     
     jjj=jjj+1;
    out_x=xk(N);
end
r1=round(r1*100/M);
r2=round(r2*100/M/L);
filter=w;
end

function [maxvalue]=MaxL2nofZeroblock(l2nwk,Zlen)
global Zi;
maxvalue=0;         
for i=1:Zlen
    index=Zi(i);
    if (maxvalue<l2nwk(index))
        maxvalue=l2nwk(index);
    end
end
end

function [Alen,Zlenout,nz,flag_NLMS]=InitialActiveBlock(M,L,thre,Zlenin,p,Alenin,nzin)
global Zi;
global w;
global Ai;
flag_NLMS=0;
l2nwk=L2nofEachBlock(M,L,w);
maxvalue=max(l2nwk);
[minvalue,minindex]=min(l2nwk);

Zlenout=Zlenin;
nz=nzin;
Alen=Alenin;
if ((Zlenin<1)||(Alenin<1))

        Zlenout=1;
        Zi(Zlenout)=minindex;  
end
if (minvalue>maxvalue/2)
     flag_NLMS=100;   
else
    [Alen,Zlenout]=LocateActiveBlock(Zlenout,L,M);    
    nz=Calculate_nz(p,Zlenout); 
    if (Alen>thre)
        flag_NLMS=100;
    else         
        v=MaxL2nofZeroblock(l2nwk,Zlenout);
        if (v>maxvalue/2)
            flag_NLMS=100; 
        else
            flag_NLMS=0;
        end
    end
end
end  

function [m]=CalculateMisalignment(nu,h,L,den)
m=0;
global w;
global Ai;
temp=zeros(L,1);
for i=1:nu
    blockindex=Ai(i);
    findex=(blockindex-1)*L+1;
    lindex= blockindex*L;
    temp=w(findex:lindex)-h(findex:lindex);
    m=m+sum(temp.*temp);
end
m=10*log10(m/den);
end


function []=UpdateCoefficients(mu,ek,v2,xk,nu,L)
temp=mu*ek/v2;
global w;
global Ai;
y=0;
for i=1:nu
    blockindex=Ai(i);
    findex=(blockindex-1)*L+1;
    lindex= blockindex*L;
    w(findex:lindex)=w(findex:lindex)+temp*xk(findex:lindex);
end
end

function   [v2]=GenerateDenominatorSPU(nu)
v2=0;
global l2n;
global Ai;
for i=1:nu
    blockindex=Ai(i);
    v2=v2+l2n(blockindex);
end
end

function [y]=GenerateOutput(nu,xk,L)
global w;
global Ai;
y=0;
for i=1:nu
    blockindex=Ai(i);
    findex=(blockindex-1)*L+1;
    lindex= blockindex*L;
    y=y+sum(xk(findex:lindex).*w(findex:lindex));    
end
end

function []=maintainL2n(M,L,out_x,xk)
global l2n;
l2n(1)=l2n(1)+xk(1)*xk(1);
for i=1:M-1
    temp=xk(i*L+1);
    temp=temp*temp;
    l2n(i)=l2n(i)-temp;
    l2n(i+1)=l2n(i+1)+temp;      
end
l2n(M)=l2n(M)-out_x*out_x;
end

function [nu]=GenerateAdaptingBlocksIndex(j,Alen,Zlen,nz)
global Ai;
global Zi;
if ((j+nz-1)>Zlen)
    nz=Zlen-j+1;
end
Ai(Alen+1:Alen+nz)=Zi(j:j+nz-1);
nu=Alen+nz;
end

function [nz]=Calculate_nz(p,Zlen)
if (p>0)
    nz=ceil(p*Zlen);
else
    if (Zlen==0)
        nz=0;
    else
        nz=1;
    end
end
end


function  [Alen,Zlen,n]=LocateActiveBlock(Zlen,L,M)
global w;
global Ai;
global Zi;
var=0;
n=0;
for i=1:Zlen
    blockindex=Zi(i);
    findex=(blockindex-1)*L+1;
    lindex= blockindex*L;
    var=var+sum(w(findex:lindex).*w(findex:lindex));
end
var=3*sqrt(var/Zlen/L);

Alen=0;
Zlen=0;
for i=1:M
    flag=0;
    findex=(i-1)*L+1;
    lindex= i*L;
    for j=findex:lindex        
        if (abs(w(j))>var)
            flag=1;
            break;
        end
        n=n+1;
    end
    if (flag==1)
        Alen=Alen+1;
        Ai(Alen)=i;
    else
        Zlen=Zlen+1;
        Zi(Zlen)=i;
    end    
end
end

function [vec]=L2nofEachBlock(M,L,xk)
vec=zeros(1,M);
for i=1:M
    findex=(i-1)*L+1;%The firtst index of the selected block;
    lindex=i*L;%The last index of the selected block; 
    vec(i)=sum(xk(findex:lindex).*xk(findex:lindex));
end
end