 function [Ynorm_actual,Ynorm_Cexact,Ynorm_Cest, e_Cexact,e_Cest, C_norm,Cest_norm, e_Cmatrix]=modelData_pest(T,n,q,p)
rng(7);

fprintf("Day %.0f\n",T-1);


V=0.5*randn(n);
Ve=V*V'; %Symmetric matrix
Ve=Ve+Ve'; %Positive definite

k=q*(p+1);
Iq=eye(q);
Ik=eye(k);

A=[zeros(q,k-q),zeros(q,q); eye(q*p),zeros(k-q,q)];
%check size of A -> kxk

B=eye(k,q);% B is a kxq matrix

Var_Uo=Ik;
Var_U(:,:,1)=Var_Uo; %A(:,:,1) means: all rows and all columns of A that are in its first page
E=zeros(size(V*randn(n,1),1),T-1);    
u=zeros(size(Iq*randn(q,1),1)   ,T-1);

for j=1:T-1-1
    E(:,j)=V*randn(n,1);    %matrix size -> nx1 for times 1 up to T-1
    u(:,j)=Iq*randn(q,1);   %matrix size -> qx1 for times 1 up to T-1 
end
u(:,T-1)=Iq*randn(q,1);       %matrix size -> qx1 at time T 

Xo=[u(:,1); zeros(q*p,1)];  %matrix size -> kx1 
X(:,1)=Xo;
Xhat0=0*Xo;                 %initial Xhat0 size kx1

Xhat(:,1)=Xhat0;            %initialise Xhat(:,1)
Xhat_t(:,1)=Xhat0;          %initialise Xhat_t(:,1)
Y=zeros(size(E,1),T-1);       %initialise Y returns -> nx1 for T periods 
C=randn(n,k);
C_norm=norm(C);
Yactual=[];
Yhat_Cexact=[];


for t=2:T%periods 
    
    %Note that t=t+1 here i.e. t is the next period and t-1 is the current
    X(:,t)=A*X(:,t-1) + B*u(:,t-1);  %computing factors of current period 
    Y(:,t-1)=C*X(:,t-1) + E(:,t-1);  %computing returns of current periods
    
    L=Var_U(:,:,t-1)*C'*inv(C*Var_U(:,:,t-1)*C'+Ve);
    
    Var_U(:,:,t)=A*Var_U(:,:,t-1)*A' +B*Iq*B' - A*L*C*Var_U(:,:,t-1)*A';
    
    %Kalman Filter
    Xhat(:,t)=A*L*(Y(:,t-1)-C*Xhat(:,t-1))+A*Xhat(:,t-1);%factors estimate
    Xhat_t(:,t-1)=Xhat(:,t-1)+L*(Y(:,t-1)-C*Xhat(:,t-1));
    
   
    if t==T
        Yactual=C*X(:,end-1); %real value of Y for next time period
        Yhat_Cexact=C*Xhat_t(:,end);%predicted value of Y using randn B
    end
end
size(X)
Xhat_Cest(:,1)=Xhat0; %initialise Xhat estimate initial
Xhat_t_Cest(:,1)=Xhat0;
Var_U_Cest(:,:,1)=Ik; %initialise covariance matric for Cest


%stacking U factors and Y returns
Uu=[];
for l=1:T-1
    Uu=[Uu;u(:,l)];
end

Yy=[];
for l=1:T-1
    Yy=[Yy;Y(:,l)];
end

size(Uu)
size(Yy)
%estimating C
Cest=Cestimate([u(:,1); zeros(q*p,1)],Uu,Yy,T-1,n,k,A,B);
Cest_norm=norm(Cest(:,end));

%adding uncertainty 
delta=0; %perturbation constant id delta=0 no perturbation injected
Cest= Cest + delta*randn(size(Cest));

Yhat_Cest=[];
for t=2:T %periods 
   
    %L using Cestimate
    L_Cest=Var_U_Cest(:,:,t-1)*Cest'*inv(Cest*Var_U_Cest(:,:,t-1)*Cest'+Ve);
    
    %Covariance using Cestimate
    Var_U_Cest(:,:,t)=A*Var_U_Cest(:,:,t-1)*A'+B*Iq*B' - A*L_Cest*Cest*Var_U_Cest(:,:,t-1)*A';
    
    %Kalman Filter using Cestimate
    Xhat_Cest(:,t)=A*L_Cest*(Y(:,t-1)-Cest*Xhat_Cest(:,t-1))+A*Xhat_Cest(:,t-1);
    Xhat_t_Cest(:,t-1)=Xhat_Cest(:,t-1)+L_Cest*(Y(:,t-1)-Cest*Xhat_Cest(:,t-1));
    
    if t==T
        Yhat_Cest=Cest*Xhat_t_Cest(:,t-1); %
        %Ynext_period=C*X(:,t); %Y next period with known C
     end 
end

Ynorm_actual=norm(Yactual(:,end)); %real value of Y for next time period
Ynorm_Cexact=norm(Yhat_Cexact(:,end)); %predicted value of Y using randn B
Ynorm_Cest=norm(Yhat_Cest);

e_Cexact=norm(Yactual(:,end)-Yhat_Cexact(:,end))/norm(Yactual(:,end));
e_Cest=norm(Yactual(:,end)-Yhat_Cest)/norm(Yactual(:,end));

e_Cmatrix=norm(C(:,end)-Cest(:,end))/norm(C(:,end));
end