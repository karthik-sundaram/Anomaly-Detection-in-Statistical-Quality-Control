d1=0;d2=0;d3=0;d4=0;d5=0;d6=0;d7=0;d8=0;d9=0;d10=0;d11=0;d12=0;d13=0;d14=0;d15=0;d16=0;d17=0;d18=0;d19=0;d20=0;d21=0;d22=0;d23=0;d24=0;d25=0;
for i=1:25
    Projectdata(1001,i)=0;
    for j=1:1000
        Projectdata(1001,i)=Projectdata(1001,i)+Projectdata(j,i);
    end
     Projectdata(1001,i)=Projectdata(1001,i)/1000;
end
for i=1:25
    for j=1:1000
        B(j,i)=Projectdata(j,i)-Projectdata(1001,i); %% B is a 1000*25 matrix
    end
end
D=transpose(B);
c=zeros(25,25);
x=zeros(25,25);
for i=1:1000
    x=D(:,i)*B(i,:);
    c=c+x;
end
c1=c/999;
s=trace(c1);
[d,v]=eig(c1);
f=trace(v);
m=1000;
n=25;
d=d(:,n:-1:1);
v=v(n:-1:1,n:-1:1);
lambda=diag(v);
%%for i=1:1000
%%t(i,1)=B(i,:)*c1*D(:,i);
    %t(i,1)=B(i,:)*D(:,i)*B(i,:)*D(:,i);
%%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mdl=zeros(1,100);
ai=zeros(1,100); 
gi=zeros(1,100);
for i=1:1:n 
    gi(i)=geomean(lambda(i:n,1)); 
    ai(i)=mean(lambda(i:n,1)); 
end; 
for j=1:1:n 
    mdl(j)=(m*(n-(j-1))*log(ai(j)/gi(j)))+(((j-1)*((2*n)-(j-1))*log(m))/2);
end;
k=0;
figure 
plot(mdl,'-OK')
hold on 
xL=get(gca,'XLim'); 
xlabel('data index')
ylabel('MDL values') 
hold off;
for j=1:1:n-1 
    if mdl(j+1)<mdl(j)
        k=j+1;
    else k; 
    end; 
end;
figure
plot(lambda,'-OK')
hold on 
xlabel('data index')
ylabel('eigen values')
hold off;

for i= 1:25
    for j=1:1000
        A(j,i)=Projectdata(j,i);
    end
end
mdata=A'; 
mu=mean(mdata,2);
y1=d(:,1)'*(mdata-repmat(mu,1,m)); 
y2=d(:,2)'*(mdata-repmat(mu,1,m));
y3=d(:,3)'*(mdata-repmat(mu,1,m));
y4=d(:,4)'*(mdata-repmat(mu,1,m)); 
y5=d(:,5)'*(mdata-repmat(mu,1,m)); 
y6=d(:,6)'*(mdata-repmat(mu,1,m));
y7=d(:,7)'*(mdata-repmat(mu,1,m));
y8=d(:,8)'*(mdata-repmat(mu,1,m)); 
y9=d(:,9)'*(mdata-repmat(mu,1,m)); 
y10=d(:,10)'*(mdata-repmat(mu,1,m));
y11=d(:,11)'*(mdata-repmat(mu,1,m));
y12=d(:,12)'*(mdata-repmat(mu,1,m)); 
y13=d(:,13)'*(mdata-repmat(mu,1,m)); 
y14=d(:,14)'*(mdata-repmat(mu,1,m));
y15=d(:,15)'*(mdata-repmat(mu,1,m));
y16=d(:,16)'*(mdata-repmat(mu,1,m)); 
y17=d(:,17)'*(mdata-repmat(mu,1,m)); 
y18=d(:,18)'*(mdata-repmat(mu,1,m));
y19=d(:,19)'*(mdata-repmat(mu,1,m));
y20=d(:,20)'*(mdata-repmat(mu,1,m)); 
y21=d(:,21)'*(mdata-repmat(mu,1,m)); 
y22=d(:,22)'*(mdata-repmat(mu,1,m));
y23=d(:,23)'*(mdata-repmat(mu,1,m));
y24=d(:,24)'*(mdata-repmat(mu,1,m)); 
y25=d(:,25)'*(mdata-repmat(mu,1,m));

UCL1=3.85*sqrt(var(y1));
LCL1=-3.85*sqrt(var(y1));
UCL2=3.85*sqrt(var(y2));
LCL2=-3.85*sqrt(var(y2));
UCL3=3.85*sqrt(var(y3)); 
LCL3=-3.85*sqrt(var(y3));
UCL4=3.85*sqrt(var(y4)); 
LCL4=-3.85*sqrt(var(y4));
UCL5=3.85*sqrt(var(y5));
LCL5=-3.85*sqrt(var(y5));
UCL6=3.85*sqrt(var(y6));
LCL6=-3.85*sqrt(var(y6));
UCL7=3.85*sqrt(var(y7)); 
LCL7=-3.85*sqrt(var(y7));
UCL8=3.85*sqrt(var(y8)); 
LCL8=-3.85*sqrt(var(y8));
UCL9=3.85*sqrt(var(y9));
LCL9=-3.85*sqrt(var(y9));
UCL10=3.85*sqrt(var(y10));
LCL10=-3.85*sqrt(var(y10));
UCL11=3.85*sqrt(var(y11)); 
LCL11=-3.85*sqrt(var(y11));
UCL12=3.85*sqrt(var(y12)); 
LCL12=-3.85*sqrt(var(y12));
UCL13=3.85*sqrt(var(y13));
LCL13=-3.85*sqrt(var(y13));
UCL14=3.85*sqrt(var(y14));
LCL14=-3.85*sqrt(var(y14));
UCL15=3.85*sqrt(var(y15)); 
LCL15=-3.85*sqrt(var(y15));
UCL16=3.85*sqrt(var(y16)); 
LCL16=-3.85*sqrt(var(y16));
UCL17=3.85*sqrt(var(y17));
LCL17=-3.85*sqrt(var(y17));
UCL18=3.85*sqrt(var(y18));
LCL18=-3.85*sqrt(var(y18));
UCL19=3.85*sqrt(var(y19)); 
LCL19=-3.85*sqrt(var(y19));
UCL20=3.85*sqrt(var(y20)); 
LCL20=-3.85*sqrt(var(y20));
UCL21=3.85*sqrt(var(y21));
LCL21=-3.85*sqrt(var(y21));
UCL22=3.85*sqrt(var(y22));
LCL22=-3.85*sqrt(var(y22));
UCL23=3.85*sqrt(var(y23)); 
LCL23=-3.85*sqrt(var(y23));
UCL24=3.85*sqrt(var(y24)); 
LCL24=-3.85*sqrt(var(y24));
UCL25=3.85*sqrt(var(y25));
LCL25=-3.85*sqrt(var(y25));

for i=1:1:m 
    if y1(i)>UCL1
        d1=d1+1;
    else if y1(i)<LCL1
            d1=d1+1;
        else d1;
        end;
    end;
end;
for i=1:1:m
    if y2(i)>UCL2
        d2=d2+1; 
    else if y2(i)<LCL2
            d2=d2+1;
else d2;
        end;
    end;
end;
for i=1:1:m
    if y3(i)>UCL3
        d3=d3+1;
    else if y3(i)<LCL3
            d3=d3+1; 
        else d3;
        end;
    end;
end;
for i=1:1:m 
    if y4(i)>UCL4
        d4=d4+1;
    else if y4(i)<LCL4
            d4=d4+1;
else d4;
end;
end; 
end; 
for i=1:1:m 
    if y5(i)>UCL5
        d5=d5+1;
    else if y5(i)<LCL5
            d5=d5+1;
        else d5;
        end;
    end;
end;
for i=1:1:m
    if y6(i)>UCL6
        d6=d6+1; 
    else if y6(i)<LCL6
            d6=d6+1;
else d6;
        end;
    end;
end;
for i=1:1:m
    if y7(i)>UCL7
        d7=d7+1;
    else if y7(i)<LCL7
            d7=d7+1; 
        else d7;
        end;
    end;
end;
for i=1:1:m 
    if y8(i)>UCL8
        d8=d8+1;
    else if y8(i)<LCL8
            d8=d8+1;
else d8;
end;
end; 
end; 
for i=1:1:m 
    if y9(i)>UCL9
        d9=d9+1;
    else if y9(i)<LCL9
            d9=d9+1;
        else d9;
        end;
    end;
end;
for i=1:1:m
    if y10(i)>UCL10
        d10=d10+1; 
    else if y10(i)<LCL10
            d10=d10+1;
else d10;
        end;
    end;
end;
for i=1:1:m
    if y11(i)>UCL11
        d11=d11+1;
    else if y11(i)<LCL11
            d11=d11+1; 
        else d11;
        end;
    end;
end;
for i=1:1:m 
    if y12(i)>UCL12
        d12=d12+1;
    else if y12(i)<LCL12
            d12=d12+1;
else d12;
end;
end; 
end; 
for i=1:1:m 
    if y13(i)>UCL13
        d13=d13+1;
    else if y13(i)<LCL13
            d13=d13+1;
        else d13;
        end;
    end;
end;
for i=1:1:m
    if y14(i)>UCL14
        d14=d14+1; 
    else if y14(i)<LCL14
            d14=d14+1;
else d14;
        end;
    end;
end;
for i=1:1:m
    if y15(i)>UCL15
        d15=d15+1;
    else if y15(i)<LCL15
            d15=d15+1; 
        else d15;
        end;
    end;
end;
for i=1:1:m 
    if y16(i)>UCL16
        d16=d16+1;
    else if y16(i)<LCL16
            d16=d16+1;
else d16;
end;
end; 
end; 
for i=1:1:m 
    if y17(i)>UCL17
        d17=d17+1;
    else if y17(i)<LCL17
            d17=d17+1;
        else d17;
        end;
    end;
end;
for i=1:1:m
    if y18(i)>UCL18
        d18=d18+1; 
    else if y18(i)<LCL18
            d18=d18+1;
else d18;
        end;
    end;
end;
for i=1:1:m
    if y19(i)>UCL19
        d19=d19+1;
    else if y19(i)<LCL19
            d19=d19+1; 
        else d19;
        end;
    end;
end;
for i=1:1:m 
    if y20(i)>UCL20
        d20=d20+1;
    else if y20(i)<LCL20
            d20=d20+1;
else d20;
end;
end; 
end; 
for i=1:1:m 
    if y21(i)>UCL21
        d21=d21+1;
    else if y21(i)<LCL21
            d21=d21+1;
        else d21;
        end;
    end;
end;
for i=1:1:m
    if y22(i)>UCL22
        d22=d22+1; 
    else if y22(i)<LCL22
            d22=d22+1;
else d22;
        end;
    end;
end;
for i=1:1:m
    if y23(i)>UCL23
        d23=d23+1;
    else if y23(i)<LCL23
            d23=d23+1; 
        else d23;
        end;
    end;
end;
for i=1:1:m 
    if y24(i)>UCL24
        d24=d24+1;
    else if y24(i)<LCL24
            d24=d24+1;
else d24;
end;
end; 
end; 
for i=1:1:m 
    if y25(i)>UCL25
        d25=d25+1;
    else if y25(i)<LCL25
            d25=d25+1;
else d25;
end;
end; 
end; 
OC1=zeros(d1,1);
k1=0; 
OC2=zeros(d2,1);
k2=0; 
OC3=zeros(d3,1); 
k3=0; 
OC4=zeros(d4,1);
k4=0;
OC5=zeros(d5,1);
k5=0; 
OC6=zeros(d6,1);
k6=0; 
OC7=zeros(d7,1); 
k7=0; 
OC8=zeros(d8,1);
k8=0;
OC9=zeros(d9,1);
k9=0; 
OC10=zeros(d10,1);
k10=0; 
OC11=zeros(d11,1); 
k11=0; 
OC12=zeros(d12,1);
k12=0;
OC13=zeros(d13,1);
k13=0; 
OC14=zeros(d14,1);
k14=0; 
OC15=zeros(d15,1); 
k15=0; 
OC16=zeros(d16,1);
k16=0;
OC17=zeros(d17,1);
k17=0; 
OC18=zeros(d18,1);
k18=0;
OC19=zeros(d19,1); 
k19=0; 
OC20=zeros(d20,1);
k20=0;
OC21=zeros(d21,1);
k21=0; 
OC22=zeros(d22,1);
k22=0; 
OC23=zeros(d23,1); 
k23=0; 
OC24=zeros(d24,1);
k24=0;
OC25=zeros(d25,1);
k25=0;


for i=1:1:m 
if y1(i)>UCL1 
k1=k1+1;
OC1(k1)=i;
else if y1(i)<LCL1 
    k1=k1+1;
    OC1(k1)=i; 
end;
end;
if y2(i)>UCL2 
k2=k2+1; 
OC2(k2)=i; 
else if y2(i)<LCL2 
    k2=k2+1; 
    OC2(k2)=i;
end;
end;
if y3(i)>UCL3 
k3=k3+1; 
OC3(k3)=i; 
else if y3(i)<LCL3 
    k3=k3+1; 
    OC3(k3)=i; 
end;
end;
if y4(i)>UCL4 
k4=k4+1;

OC4(k4)=i;
else if y4(i)<LCL4
    k4=k4+1; 
    OC4(k4)=i; 
end; 
end;
if y5(i)>UCL5 
k5=k5+1;
OC5(k5)=i;
else if y5(i)<LCL5 
    k5=k5+1;
    OC5(k5)=i; 
end;
end;
if y6(i)>UCL6 
k6=k6+1; 
OC6(k6)=i; 
else if y6(i)<LCL6 
    k6=k6+1; 
    OC6(k6)=i;
end;
end;
if y7(i)>UCL7 
k7=k7+1; 
OC7(k7)=i; 
else if y7(i)<LCL7 
    k7=k7+1; 
    OC7(k7)=i; 
end;
end;
if y8(i)>UCL8 
k8=k8+1;

OC9(k9)=i;
else if y9(i)<LCL9
    k9=k9+1; 
    OC9(k9)=i; 
end; 
end;
if y10(i)>UCL10 
k10=k10+1;
OC10(k10)=i;
else if y10(i)<LCL10 
    k10=k10+1;
    OC10(k10)=i; 
end;
end;
if y11(i)>UCL11 
k11=k11+1; 
OC11(k11)=i; 
else if y11(i)<LCL11 
    k11=k11+1; 
    OC11(k11)=i;
end;
end;
if y12(i)>UCL12 
k12=k12+1; 
OC12(k12)=i; 
else if y12(i)<LCL12 
    k12=k12+1; 
    OC12(k12)=i; 
end;
end;
if y13(i)>UCL13 
k13=k13+1;

OC13(k13)=i;
else if y13(i)<LCL13
    k13=k13+1; 
    OC13(k13)=i; 
end; 
end;
if y14(i)>UCL14 
k14=k14+1;
OC14(k14)=i;
else if y14(i)<LCL14 
    k14=k14+1;
    OC14(k14)=i; 
end;
end;
if y15(i)>UCL15 
k15=k15+1; 
OC15(k15)=i; 
else if y15(i)<LCL15 
    k15=k15+1; 
    OC15(k15)=i;
end;
end;
  if y16(i)>UCL16
        k16=k16+1; 
        OC16(k16)=i;
    else
        if y16(i)<LCL16
            k16=k16+1; 
            OC16(k16)=i;
        end 
    end
  
     if y17(i)>UCL17
        k17=k17+1; 
        OC17(k17)=i;
    else
        if y17(i)<LCL17
            k17=k17+1; 
            OC17(k17)=i;
        end 
     end
    
      if y18(i)>UCL18
        k18=k18+1; 
        OC18(k18)=i;
    else
        if y18(i)<LCL18
            k18=k18+1; 
            OC18(k18)=i;
        end 
      end
    
       if y19(i)>UCL19
        k19=k19+1; 
        OC19(k19)=i;
    else
        if y19(i)<LCL19
            k19=k19+1; 
            OC19(k19)=i;
        end 
       end
    
        if y20(i)>UCL20
        k20=k20+1; 
        OC20(k20)=i;
    else
        if y20(i)<LCL20
            k20=k20+1; 
            OC20(k20)=i;
        end 
        end
    
          if y21(i)>UCL21
        k21=k21+1; 
        OC2(k21)=i;
    else
        if y21(i)<LCL21
            k21=k21+1; 
            OC21(k21)=i;
        end 
          end
      if y22(i)>UCL22
        k22=k22+1; 
        OC22(k22)=i;
    else
        if y22(i)<LCL22
            k22=k22+1; 
            OC22(k22)=i;
        end 
      end
      if y23(i)>UCL23
        k23=k23+1; 
        OC23(k23)=i;
    else
        if y23(i)<LCL23
            k23=k23+1; 
            OC2(k23)=i;
        end 
      end
      if y24(i)>UCL24
        k24=k24+1; 
        OC24(k2)=i;
    else
        if y24(i)<LCL24
            k24=k24+1; 
            OC24(k24)=i;
        end 
      end
      if y25(i)>UCL25
        k25=k25+1; 
        OC25(k25)=i;
    else
        if y25(i)<LCL25
            k25=k25+1; 
            OC25(k25)=i;
        end 
    end

end; 
figure
plot(y1,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL1 UCL1])
line(xL,[LCL1 LCL1]) 
text(m+.5,UCL1,'UCL1') 
text(m+.5,LCL1,'LCL1') 
hold off; 
figure 
plot(y2,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL2 UCL2]) 
line(xL,[LCL2 LCL2])
text(m+.5,UCL2,'UCL2')
text(m+.5,LCL2,'LCL2')
hold off;
figure
plot(y3,'-OK') 
hold on xL=get(gca,'XLim'); 
line(xL,[UCL3 UCL3]) 
line(xL,[LCL3 LCL3]) 
text(m+.5,UCL3,'UCL3') 
text(m+.5,LCL3,'LCL3') 
hold off; 
figure 
plot(y4,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL4 UCL4])
line(xL,[LCL4 LCL4]) 
text(m+.5,UCL4,'UCL4')
text(m+.5,LCL4,'LCL4')
hold off;
figure
plot(y5,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL5 UCL5])
line(xL,[LCL5 LCL5]) 
text(m+.5,UCL5,'UCL5') 
text(m+.5,LCL5,'LCL5') 
hold off; 
figure 
plot(y6,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL6 UCL6]) 
line(xL,[LCL6 LCL6])
text(m+.5,UCL6,'UCL6')
text(m+.5,LCL6,'LCL6')
hold off;
figure 
plot(y7,'-OK') 
hold on xL=get(gca,'XLim'); 
line(xL,[UCL7 UCL7]) 
line(xL,[LCL7 LCL7]) 
text(m+.5,UCL7,'UCL7') 
text(m+.5,LCL7,'LCL7') 
hold off; 
figure
plot(y8,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL8 UCL8])
line(xL,[LCL8 LCL8]) 
text(m+.5,UCL8,'UCL8')
text(m+.5,LCL8,'LCL8')
hold off;
figure
plot(y9,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL9 UCL9])
line(xL,[LCL9 LCL9]) 
text(m+.5,UCL9,'UCL9') 
text(m+.5,LCL9,'LCL9') 
hold off; 
figure 
plot(y10,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL10 UCL10]) 
line(xL,[LCL10 LCL10])
text(m+.5,UCL10,'UCL10')
text(m+.5,LCL10,'LCL10')
hold off;
figure 
plot(y11,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL11 UCL11]) 
line(xL,[LCL11 LCL11]) 
text(m+.5,UCL11,'UCL11') 
text(m+.5,LCL11,'LCL11') 
hold off; 
figure
plot(y12,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL12 UCL12])
line(xL,[LCL12 LCL12]) 
text(m+.5,UCL12,'UCL12')
text(m+.5,LCL12,'LCL12')
hold off;
figure 
plot(y13,'-OK') 
hold on 
xL=get(gca,'XLim'); 
line(xL,[UCL13 UCL13])
line(xL,[LCL13 LCL13]) 
text(m+.5,UCL13,'UCL13') 
text(m+.5,LCL13,'LCL13') 
hold off; 
figure
plot(y14,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL14 UCL14])
line(xL,[LCL14 LCL14])
text(m+.5,UCL14,'UCL14') 
text(m+.5,LCL14,'LCL14') 
hold off;

figure
plot(y15,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL15 UCL15])
line(xL,[LCL15 LCL15])
text(m+.5,UCL15,'UCL15') 
text(m+.5,LCL15,'LCL15') 
hold off;

figure
plot(y16,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL16 UCL16])
line(xL,[LCL16 LCL16])
text(m+.5,UCL16,'UCL16') 
text(m+.5,LCL16,'LCL16') 
hold off;

figure
plot(y17,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL17 UCL17])
line(xL,[LCL17 LCL17])
text(m+.5,UCL17,'UCL17') 
text(m+.5,LCL17,'LCL17') 
hold off;

figure
plot(y18,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL18 UCL18])
line(xL,[LCL18 LCL18])
text(m+.5,UCL18,'UCL18') 
text(m+.5,LCL18,'LCL18') 
hold off;

figure
plot(y19,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL19 UCL19])
line(xL,[LCL19 LCL19])
text(m+.5,UCL19,'UCL19') 
text(m+.5,LCL19,'LCL19') 
hold off;

figure
plot(y20,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL20 UCL20])
line(xL,[LCL20 LCL20])
text(m+.5,UCL20,'UCL20') 
text(m+.5,LCL20,'LCL20') 
hold off;


figure
plot(y21,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL21 UCL21])
line(xL,[LCL21 LCL21])
text(m+.5,UCL21,'UCL21') 
text(m+.5,LCL21,'LCL21') 
hold off;

figure
plot(y22,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL22 UCL22])
line(xL,[LCL22 LCL22])
text(m+.5,UCL22,'UCL22') 
text(m+.5,LCL22,'LCL22') 
hold off;

figure
plot(y23,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL23 UCL23])
line(xL,[LCL23 LCL23])
text(m+.5,UCL23,'UCL23') 
text(m+.5,LCL23,'LCL23') 
hold off;


figure
plot(y24,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL24 UCL24])
line(xL,[LCL24 LCL24])
text(m+.5,UCL24,'UCL24') 
text(m+.5,LCL24,'LCL24') 
hold off;

figure
plot(y25,'-OK')
hold on
xL=get(gca,'XLim');
line(xL,[UCL25 UCL25])
line(xL,[LCL25 LCL25])
text(m+.5,UCL25,'UCL25') 
text(m+.5,LCL25,'LCL25') 
hold off;

v1=0;
for i=1:1:25
    v1=v1+v(i,i);
end
for i=1:1:25
    vp(i,i)=v(i,i)/v1*100;
end
for i=1:1:25
    vp1(i,i)=0;
    for j=1:1:i
        vp1(i,i)=vp1(i,i)+vp(j,j);
    end
    end

Projectdata1=Projectdata([1:370,372:766,768:790,792:797,799:823,825:838,840:842,844:846,848:878,880:897,899:952,954:end],:)
