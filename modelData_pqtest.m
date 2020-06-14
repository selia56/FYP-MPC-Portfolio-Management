function [Ynorm_actual,Ynorm_Cexact,Ynorm_Cest,e_Cexact,e_norm_Cest,epsilon]=modelData_pqtest(T,n,q_actual,p_actual,N,q_test,p_test,q_est)

%% 
rng(7);

V=0.2*randn(n);
Ve=V*V'; %Symmetric matrix
Ve=Ve+Ve'; %Positive definite

k=q_actual*(p_actual+1);
Iq=eye(q_actual);
Ik=eye(k);

A=[zeros(q_actual,k-q_actual),zeros(q_actual,q_actual); eye(q_actual*p_actual),zeros(k-q_actual,q_actual)];%kxk

B=eye(k,q_actual);% B is a kxq matrix

Var_Uo=Ik;
Var_U(:,:,1)=Var_Uo; %A(:,:,1) means: all rows and all columns of A that are in its first page kxk
E=zeros(size(V*randn(n,1),1),T-1); 

u=zeros(size(Iq*randn(q_actual,1),1),T-1);


for j=1:T-1-1
   
    E(:,j)=V*randn(n,1);    %matrix size -> nx1 for times 1 up to T-1
    
    u(:,j)=Iq*randn(q_actual,1);   %matrix size -> qx1 for times 1 up to T-1 
end


u(:,T-1)=Iq*randn(q_actual,1);       %matrix size -> qx1 at time T 
Xo=[u(:,1); zeros(q_actual*p_actual,1)];  %matrix size -> kx1 
X(:,1)=Xo;

Xhat0=0*Xo;                 %initial Xhat0 size kx1

Xhat(:,1)=Xhat0;            %initialise Xhat(:,1)
Xhat_t(:,1)=Xhat0;          %initialise Xhat_t(:,1)
Y=zeros(size(E,1),T-1);       %initialise Y returns -> nx1 for T periods 
C=randn(n,k);
C_norm=norm(C);

Yactual=[];
Yhat_Cq_exact=[];

fprintf("Day %.0f\n",T-1);
%Kalman Filter tracking Factor Matrix X_{t+1|t}
for t=2:T%periods 
    
    %Note that t=t+1 here i.e. t is the next period and t-1 is the current
    X(:,t)=A*X(:,t-1) + B*u(:,t-1);  %computing factors of current period 
    Y(:,t-1)=C*X(:,t-1) + E(:,t-1);  %computing returns of current periods
    
    L_q=Var_U(:,:,t-1)*C'*inv(C*Var_U(:,:,t-1)*C'+Ve);
    
    Var_U(:,:,t)=A*Var_U(:,:,t-1)*A' +B*Iq*B' - A*L_q*C*Var_U(:,:,t-1)*A';
    
    %Kalman Filter
    Xhat(:,t)=A*L_q*(Y(:,t-1)-C*Xhat(:,t-1))+A*Xhat(:,t-1); %factors estimate
    Xhat_t(:,t-1)=Xhat(:,t-1)+L_q*(Y(:,t-1)-C*Xhat(:,t-1));
    
   
    if t==T
        Yactual=C*X(:,end-1); %real value of Y for next time period
        Yhat_Cq_exact=C*Xhat_t(:,end);%predicted value of Y using exact C
    end
end
%%

%q-estimate

Yy=[];
    for l=1:T-1
        Yy=[Yy;Y(:,l)];
    end
Yhat_Cq_est=zeros(n,q_test);

%%


%q_est=1;    
 for p_est=1:p_test
Uu=[]; 
  Xhat_q_est=[];
  Xhat_t_Cq_est=[];
  Var_U_q_est=[];
  
    
    k_q_est=q_est*(p_est+1); %initialise k_est
   
    Iq_est=eye(q_est);
    Iq_test=eye(q_test);
    
    
    u_q=zeros(size(Iq_test*randn(q_test,1),1),T-1);
    u_qq=zeros(size(Iq_test*randn(q_test,1),1),T-1);
    
    for j=1:T-1-1
        u_q(:,j)=Iq_test*randn(q_test,1);   %matrix size -> qx1 for times 1 up to T-1 
    end
    
    
    u_q(:,T-1)=Iq_test*randn(q_test,1);       %matrix size -> qx1 at time T  
    
     for j=1:T-1
        u_qq(1:size(u,1),j)=u(:,j);   %matrix size -> qx1 for times 1 up to T-1 
     end
 
    
    if size(u,1)==q_test
        break
    else 
        
      u_qq(size(u,1)+1:end,:)=eye(q_test-(size(u,1)))*randn(q_test-(size(u,1)),size(u,2));
    end
    
    X0_q_est=[u_qq(1:q_est,1); zeros(q_est*p_est,1)];
    
    Xhat_q_est(:,1)=0*X0_q_est; %initialise Xhat estimate initial
    
    A_q_est=[zeros(q_est,k_q_est-q_est),zeros(q_est,q_est); eye(q_est*p_est),zeros(k_q_est-q_est,q_est)]; %kxk
    B_q_est=eye(k_q_est,q_est);% B is a kxq matrix
   
    Xhat_t_Cq_est(:,1)=0*X0_q_est;
   
    Var_U_q_est(:,:,1)=eye(k_q_est); %initialise covariance matric for Cest
    Cq_est=zeros(n,q_est*(p_est+1));
    
    for l=1:T-1
        Uu=[Uu;u_qq(1:q_est,l)];
    end

    Cq_est=Cestimate([u_qq(1:q_est,1); zeros(q_est*p_est,1)],Uu,Yy,T-1,n,k_q_est,A_q_est,B_q_est);
    
   
    %fprintf("P estimate = %.0f\n",p_est);
    
    for t=2:T %periods 
    
        %L using Cestimate
        L_Cq_est=Var_U_q_est(:,:,t-1)*Cq_est'*inv(Cq_est*Var_U_q_est(:,:,t-1)*Cq_est'+Ve);

        %Covariance using Cestimate
        Var_U_q_est(:,:,t)=A_q_est*Var_U_q_est(:,:,t-1)*A_q_est'+B_q_est*Iq_est*B_q_est' - A_q_est*L_Cq_est*Cq_est*Var_U_q_est(:,:,t-1)*A_q_est';
        
        %Kalman Filter using Cestimate
        Xhat_q_est(:,t)=A_q_est*L_Cq_est*(Y(:,t-1)-Cq_est*Xhat_q_est(:,t-1))+A_q_est*Xhat_q_est(:,t-1);
        Xhat_t_Cq_est(:,t-1)=Xhat_q_est(:,t-1)+L_Cq_est*(Y(:,t-1)-Cq_est*Xhat_q_est(:,t-1));
      
        if t==T
            Yhat_Cq_est(:,p_est,:)= Cq_est*Xhat_t_Cq_est(:,t-1); %
            %Ynext_period=C*X(:,t); %Y next period with known C
        end 
    end  
    
    
    Ynorm_Cest(p_est,:,:)=norm(Yhat_Cq_est(:,p_est,:));
    
    e_norm_Cest(p_est,:,:)=norm(Yactual(:,end)-Yhat_Cq_est(:,p_est,:))/norm(Yactual(:,end));
    
    e_Cexact=norm(Yactual(:,end)-Yhat_Cq_exact(:,end))/norm(Yactual(:,end));
    Ynorm_actual=norm(Yactual); %real value of Y for next time period
    Ynorm_Cexact=norm(Yhat_Cq_exact);  
     
    epsilon=norm(Yactual(:,end)-Yhat_Cq_est(:,p_est,:))
 end


end