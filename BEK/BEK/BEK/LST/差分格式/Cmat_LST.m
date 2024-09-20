function [ A,B,C ] = Cmat_LST( NY,omge,beta,Ma,Re0,Pr,gamma,ix,kind )
load( './baseflow.mat' );
[T_L,A_L,B_L,C_L,D_L,Vxx,Vyy,Vxy,Vxz,Vyz,Vzz]=deal( zeros(5*NY+5,5*NY+5) );
[ LY1,LY2 ] = Dy( NY,2 );  [ L1,~ ] = Dy( NY,1 );
for k=1:NY+1
    m=(k-1)*5;
    T_L(m+1,m+1)=1.0;  T_L(m+2,m+2)=den(k,ix);  T_L(m+3,m+3)=den(k,ix);  T_L(m+4,m+4)=den(k,ix);
    T_L(m+5,m+1)=(1.0-gamma)/gamma*T(k,ix);  T_L(m+5,m+5)=den(k,ix)/gamma;
    A_L(m+1,m+1)=u(k,ix);  A_L(m+1,m+2)=den(k,ix);
    A_L(m+2,m+1)=T(k,ix)/gamma/Ma/Ma;  A_L(m+2,m+2)=den(k,ix)*u(k,ix);  A_L(m+2,m+3)=-dmdy(k,ix)/Re0;  A_L(m+2,m+5)=den(k,ix)/gamma/Ma/Ma;
    A_L(m+3,m+2)=2.0/3.0*dmdy(k,ix)/Re0;  A_L(m+3,m+3)=den(k,ix)*u(k,ix);  A_L(m+3,m+5)=-dmdT(k,ix)/Re0*(dudy(k,ix));
    A_L(m+4,m+4)=den(k,ix)*u(k,ix);
    A_L(m+5,m+1)=(1.0-gamma)/gamma*T(k,ix)*u(k,ix);  A_L(m+5,m+3)=-(gamma-1.0)*Ma*Ma/Re0*mul(k,ix)*2.0*(dudy(k,ix));
    A_L(m+5,m+5)=den(k,ix)*u(k,ix)/gamma;
    B_L(m+1,m+3)=den(k,ix);  B_L(m+2,m+2)=(-dmdy(k,ix)/Re0);  B_L(m+2,m+5)=-dmdT(k,ix)/Re0*dudy(k,ix);
    B_L(m+3,m+1)=T(k,ix)/gamma/Ma/Ma;  B_L(m+3,m+3)=(-4.0/3.0/Re0*dmdy(k,ix));  B_L(m+3,m+5)=den(k,ix)/Ma/Ma/gamma;
    B_L(m+4,m+4)=(-dmdy(k,ix)/Re0);  B_L(m+4,m+5)=(-dmdT(k,ix)/Re0*dwdy(k,ix));
    B_L(m+5,m+2)=(-(gamma-1.0)*Ma*Ma/Re0*2.0*mul(k,ix)*(dudy(k,ix)));
    B_L(m+5,m+4)=(-(gamma-1.0)*Ma*Ma/Re0*2.0*mul(k,ix)*(dwdy(k,ix)));
    B_L(m+5,m+5)=(-(dmdy(k,ix)+dmdT(k,ix)*dTdy(k,ix))/Re0/Pr);
    C_L(m+1,m+1)=w(k,ix);  C_L(m+1,m+4)=den(k,ix);
    C_L(m+2,m+2)=den(k,ix)*w(k,ix);
    C_L(m+3,m+3)=den(k,ix)*w(k,ix);  C_L(m+3,m+4)=2.0/3.0*dmdy(k,ix)/Re0;  C_L(m+3,m+5)=-dmdT(k,ix)/Re0*dwdy(k,ix);
    C_L(m+4,m+1)=T(k,ix)/gamma/Ma/Ma;  C_L(m+4,m+3)=-dmdy(k,ix)/Re0;  C_L(m+4,m+4)=den(k,ix)*w(k,ix);  C_L(m+4,m+5)=den(k,ix)/gamma/Ma/Ma;
    C_L(m+5,m+1)=(1-gamma)/gamma*T(k,ix)*w(k,ix);  C_L(m+5,m+3)=-(gamma-1.0)*Ma*Ma/Re0*2.0*mul(k,ix)*(dwdy(k,ix));
    C_L(m+5,m+5)=den(k,ix)*w(k,ix)/gamma;
    D_L(m+1,m+3)=dendy(k,ix);
    D_L(m+2,m+3)=den(k,ix)*dudy(k,ix);
    D_L(m+2,m+5)=-dmdTdy(k,ix)/Re0*(dudy(k,ix))-dmdT(k,ix)/Re0*(dudy2(k,ix));
    D_L(m+3,m+1)=dTdy(k,ix)/gamma/Ma/Ma;  D_L(m+3,m+5)=dendy(k,ix)/gamma/Ma/Ma;
    D_L(m+4,m+3)=den(k,ix)*dwdy(k,ix);  D_L(m+4,m+5)=-dmdTdy(k,ix)/Re0*dwdy(k,ix)-dmdT(k,ix)/Re0*dwdy2(k,ix);
    D_L(m+5,m+3)=dTdy(k,ix)*den(k,ix)/gamma+(1.0-gamma)/gamma*T(k,ix)*dendy(k,ix);
    D_L(m+5,m+5)=-1.0/Re0/Pr*(dmdTdy(k,ix)*dTdy(k,ix)+dmdT(k,ix)*dTdy2(k,ix)) ...
        -(gamma-1.0)*Ma*Ma/Re0*dmdT(k,ix)*(dudy(k,ix)*dudy(k,ix)+dwdy(k,ix)*dwdy(k,ix));
    Vxx(m+2,m+2)=-4.0/3.0*mul(k,ix)/Re0;   Vxx(m+3,m+3)=-mul(k,ix)/Re0;  Vxx(m+4,m+4)=-mul(k,ix)/Re0;  Vxx(m+5,m+5)=-mul(k,ix)/Re0/Pr;
    Vyy(m+2,m+2)=-mul(k,ix)/Re0;  Vyy(m+3,m+3)=-4.0/3.0*mul(k,ix)/Re0;  Vyy(m+4,m+4)=-mul(k,ix)/Re0;  Vyy(m+5,m+5)=-mul(k,ix)/Re0/Pr;
    Vzz(m+2,m+2)=-mul(k,ix)/Re0;   Vzz(m+3,m+3)=-mul(k,ix)/Re0;  Vzz(m+4,m+4)=-4.0/3.0*mul(k,ix)/Re0;  Vzz(m+5,m+5)=-mul(k,ix)/Re0/Pr;
    Vxy(m+2,m+3)=-1.0/3.0*mul(k,ix)/Re0;  Vxy(m+3,m+2)=-1.0/3.0*mul(k,ix)/Re0;
    Vxz(m+2,m+4)=-1.0/3.0*mul(k,ix)/Re0;  Vxz(m+4,m+2)=-1.0/3.0*mul(k,ix)/Re0;
    Vyz(m+3,m+4)=-1.0/3.0*mul(k,ix)/Re0;  Vyz(m+4,m+3)=-1.0/3.0*mul(k,ix)/Re0;
end

if(kind==1)
    Vxx(5*NY+1:5*NY+5,:)=0;  Vyy(5*NY+1:5*NY+5,:)=0;  Vxy(5*NY+1:5*NY+5,:)=0;
    Vxz(5*NY+1:5*NY+5,:)=0;  Vyz(5*NY+1:5*NY+5,:)=0;  Vzz(5*NY+1:5*NY+5,:)=0;
    TF(1:5,1:5)=T_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5);  BF(1:5,1:5)=B_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5);  
    [V,E]=eig(TF\BF);  E(1,1)=max(0,E(1,1));  E(2,2)=max(0,E(2,2));  E(3,3)=max(0,E(3,3));  E(4,4)=max(0,E(4,4));  E(5,5)=max(0,E(5,5));
    BF=TF*(V*E/V);  B_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5)=BF(:,:);
    A=D_L+complex(0,1)*beta*C_L-complex(0,1)*omge*T_L-beta*beta*Vzz+B_L*LY1+complex(0,1)*beta*Vyz*LY1+Vyy*LY2;
    B=complex(0,1)*Vxy*LY1-beta*Vxz+complex(0,1)*A_L;
    C=-Vxx;
    A=sparse(A);  B=sparse(B);  C=sparse(C);
end

if(kind==2)
    [dB1,dB2,dB3,dC1,dC2]=deal( zeros(5*NY+5,5*NY+5) );
    [T13,T22,T25,T31,T33,T35,T44,T45,T52,T54,T55]=deal( zeros(NY+1,1) );
    for k=1:NY+1
        T13(k)=den(k,ix);
        T22(k)=(-dmdy(k,ix)/Re0);
        T25(k)=-dmdT(k,ix)/Re0*dudy(k,ix);
        T31(k)=T(k,ix)/gamma/Ma/Ma;
        T33(k)=(-4.0/3.0/Re0*dmdy(k,ix));
        T35(k)=den(k,ix)/Ma/Ma/gamma;
        T44(k)=(-dmdy(k,ix)/Re0);
        T45(k)=(-dmdT(k,ix)/Re0*dwdy(k,ix));
        T52(k)=(-(gamma-1.0)*Ma*Ma/Re0*2.0*mul(k,ix)*(dudy(k,ix)));
        T54(k)=(-(gamma-1.0)*Ma*Ma/Re0*2.0*mul(k,ix)*(dwdy(k,ix)));
        T55(k)=(-(dmdy(k,ix)+dmdT(k,ix)*dTdy(k,ix))/Re0/Pr);
    end
    T13=L1*T13; T22=L1*T22; T25=L1*T25; T31=L1*T31; T33=L1*T33; T35=L1*T35; T44=L1*T44; T45=L1*T45; T52=L1*T52; T54=L1*T54; T55=L1*T55;
    for k=1:NY+1
        m=5*(k-1);
        dB1(m+1,m+3)=T13(k);  dB1(m+2,m+2)=T22(k);  dB1(m+2,m+5)=T25(k);  dB1(m+3,m+1)=T31(k);  dB1(m+3,m+3)=T33(k);
        dB1(m+3,m+5)=T35(k);  dB1(m+4,m+4)=T44(k);  dB1(m+4,m+5)=T45(k);  dB1(m+5,m+2)=T52(k);  dB1(m+5,m+4)=T54(k);  dB1(m+5,m+5)=T55(k);
    end
    clear T13 T22 T25 T31 T33 T35 T44 T45 T52 T54 T55;
    [T34,T43]=deal( zeros(NY+1,1) );
    for k=1:NY+1
        T34(k)=-1.0/3.0*mul(k,ix)/Re0;  T43(k)=-1.0/3.0*mul(k,ix)/Re0;
    end
    T34=L1*T34;  T43=L1*T43;
    for k=1:NY+1
        m=5*(k-1);  dB2(m+3,m+4)=T34(k);  dB2(m+4,m+3)=T43(k);
    end
    clear T34 T43;
    [T23,T32]=deal( zeros(NY+1,1) );
    for k=1:NY+1
        T23(k)=-1.0/3.0*mul(k,ix)/Re0;   T32(k)=-1.0/3.0*mul(k,ix)/Re0;
    end
    T23=L1*T23;  T32=L1*T32;
    for k=1:NY+1
        m=5*(k-1);  dB3(m+2,m+3)=T23(k);  dB3(m+3,m+2)=T32(k);
    end
    clear T23 T32;
    [T22,T33,T44,T55]=deal( zeros(NY+1,1) );
    for k=1:NY+1
        T22(k)=-mul(k,ix)/Re0;  T33(k)=-4.0/3.0*mul(k,ix)/Re0;
        T44(k)=-mul(k,ix)/Re0;  T55(k)=-mul(k,ix)/Re0/Pr;
    end
    T22=L1*T22; T33=L1*T33; T44=L1*T44; T55=L1*T55;
    for k=1:NY+1
        m=5*(k-1);  dC1(m+2,m+2)=T22(k);  dC1(m+3,m+3)=T33(k);  dC1(m+4,m+4)=T44(k);  dC1(m+5,m+5)=T55(k);
    end
    T22=L1*T22;  T33=L1*T33;  T44=L1*T44;  T55=L1*T55;
    for k=1:NY+1
        m=5*(k-1);  dC2(m+2,m+2)=T22(k);  dC2(m+3,m+3)=T33(k);  dC2(m+4,m+4)=T44(k);  dC2(m+5,m+5)=T55(k);
    end
    clear T22 T33 T44 T55;

    T_L=T_L.';  A_L=A_L.';  B_L=-B_L.';  C_L=C_L.';  D_L=D_L.';  Vxx=Vxx.';  Vyy=Vyy.';  Vxy=-Vxy.';  Vxz=Vxz.';  Vyz=-Vyz.';  Vzz=Vzz.';
    dB1=dB1.';  dB2=dB2.';  dB3=dB3.';  dC1=dC1.';  dC2=dC2.';
    Vxx(5*NY+1:5*NY+5,:)=0;  Vyy(5*NY+1:5*NY+5,:)=0;  Vxy(5*NY+1:5*NY+5,:)=0;  
    Vxz(5*NY+1:5*NY+5,:)=0;  Vyz(5*NY+1:5*NY+5,:)=0;  Vzz(5*NY+1:5*NY+5,:)=0;
    TF(1:5,1:5)=T_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5);  BF(1:5,1:5)=B_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5);  
    [V,E]=eig(TF\BF);  E(1,1)=max(0,E(1,1));  E(2,2)=max(0,E(2,2));  E(3,3)=max(0,E(3,3));  E(4,4)=max(0,E(4,4));  E(5,5)=max(0,E(5,5));
    BF=TF*(V*E/V);  B_L(5*NY+1:5*NY+5,5*NY+1:5*NY+5)=BF(:,:);
    A=D_L+complex(0,1)*beta*C_L-complex(0,1)*omge*T_L-beta*beta*Vzz+B_L*LY1+complex(0,1)*beta*Vyz*LY1+Vyy*LY2 ...
        +2.0*dC1*LY1-dB1-complex(0,1)*beta*dB2+dC2;
    B=complex(0,1)*Vxy*LY1-beta*Vxz+complex(0,1)*A_L-complex(0,1)*dB3;
    C=-Vxx;
    A=sparse(A);  B=sparse(B);  C=sparse(C);
end










