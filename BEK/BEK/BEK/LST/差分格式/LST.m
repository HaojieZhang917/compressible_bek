clear all;
imax=1;  [Ma,Re0,Fw,beta,gamma,Pr,NY,iadjoint,bcwall]=flow_parameters;  omge=Fw*Re0;
o_Matrix=sparse(5*(NY+1),5*(NY+1));  one_Matrix=speye(5*(NY+1),5*(NY+1));
load('D:\POLESTAR\crossflow\baseflow.mat','x','y');
eig_init=2.1-0.1i;  eig_lst=zeros(imax,1);
for ix=1:imax
    iloc=300;
    if(iadjoint==0)
        [ A,B,C ] = Cmat_LST( NY,omge,beta,Ma,Re0,Pr,gamma,iloc,1 );
    else
        [ A,B,C ] = Cmat_LST( NY,omge,beta,Ma,Re0,Pr,gamma,iloc,2 );
    end
    L_Matrix11=o_Matrix;    L_Matrix12=one_Matrix;
    L_Matrix21=A;  L_Matrix21=full(L_Matrix21);  L_Matrix21(2:5,:)=0.0;
    L_Matrix21(2,2)=1.0;  L_Matrix21(3,3)=1.0;  L_Matrix21(4,4)=1.0;
    if(bcwall==0)
        L_Matrix21(5,5)=1.0;
    else
        [ LY1,~ ] = Dy( NY,2 );  LY1=full(LY1);  L_Matrix21(5,:)=LY1(5,:);
    end
    L_Matrix21=sparse(L_Matrix21);  L_Matrix22=B;  L_Matrix22=full(L_Matrix22);  L_Matrix22(2:5,:)=0.0;  L_Matrix22=sparse(L_Matrix22);
    R_Matrix11=one_Matrix;  R_Matrix12=o_Matrix;  R_Matrix21=o_Matrix; 
    R_Matrix22=-C;  R_Matrix22=full(R_Matrix22);  R_Matrix22(2:5,:)=0.0;  R_Matrix22=sparse(R_Matrix22);
    L_Matrix=[L_Matrix11 L_Matrix12;L_Matrix21 L_Matrix22];
    R_Matrix=[R_Matrix11 R_Matrix12;R_Matrix21 R_Matrix22];
    num=1;  opts.maxit=500;  opts.tol=1e-9;  opts.isreal=false;  sigma=eig_init;
    [Vec,Val,flag]=eigs(L_Matrix,R_Matrix,num,sigma,opts);  eig_lst(ix)=Val;  eig_init=eig_lst(ix);
    disp(['Rex=',num2str(Re0),' , ','eig','(',num2str(ix),')=',num2str(eig_lst(ix),'%16.14f')]);
    [u_r,u_i,u_m,v_r,v_i,v_m,w_r,w_i,w_m,d_r,d_i,d_m,T_r,T_i,T_m,temp]=deal(zeros(NY+1,1));
    for i=1:NY+1
        temp(i)=abs(Vec(2+5*(i-1)));
    end
    max_prof=max(temp);
    for i=1:NY+1
        d_i(i)=imag(Vec(1+5*(i-1)))/max_prof;  d_r(i)=real(Vec(1+5*(i-1)))/max_prof;  d_m(i)=sqrt(d_i(i)*d_i(i)+d_r(i)*d_r(i));
        u_i(i)=imag(Vec(2+5*(i-1)))/max_prof;  u_r(i)=real(Vec(2+5*(i-1)))/max_prof;  u_m(i)=sqrt(u_i(i)*u_i(i)+u_r(i)*u_r(i));
        v_i(i)=imag(Vec(3+5*(i-1)))/max_prof;  v_r(i)=real(Vec(3+5*(i-1)))/max_prof;  v_m(i)=sqrt(v_i(i)*v_i(i)+v_r(i)*v_r(i));
        w_i(i)=imag(Vec(4+5*(i-1)))/max_prof;  w_r(i)=real(Vec(4+5*(i-1)))/max_prof;  w_m(i)=sqrt(w_i(i)*w_i(i)+w_r(i)*w_r(i));
        T_i(i)=imag(Vec(5+5*(i-1)))/max_prof;  T_r(i)=real(Vec(5+5*(i-1)))/max_prof;  T_m(i)=sqrt(T_i(i)*T_i(i)+T_r(i)*T_r(i));
    end
end
if(iadjoint==1)
    Vec_adjoint=Vec(1:5*NY+5)/max_prof;  save('D:\POLESTAR\crossflow\Vec_adjoint.mat','Vec_adjoint');
else
    Vec_normal=Vec(1:5*NY+5)/max_prof;  save('D:\POLESTAR\crossflow\Vec_normal.mat','Vec_normal');
    save('D:\POLESTAR\crossflow\eig_lst.mat','eig_lst');
end
if(iadjoint==0)
    growth=-imag(eig_lst);  Nt=zeros(imax,1);
    for ix=2:imax
        Nt(ix)=trapz(x(1:ix),growth(1:ix));
    end
end









