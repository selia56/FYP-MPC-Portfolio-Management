function Cest=Cestimate(Xo,Uu,Yy,T,n,k,A,B)

cvx_begin
cvx_begin quiet

    variable Cc(n,k);
    
    %c0
    c0=[];
    c0=kron(eye(T),Cc);          %c0=[C 0 0; 0 C 0; 0 0 C]
    
    %A0
    A0=[];
    for i=1:T
        A0=[A0; repmat(A^(i-1),1,1)];
    end
    
    %Bb
    Bb=kron(eye(T),B);
    
    %Aa
    Aa = zeros((T)*k);
   for ii = 1:k:length(Aa)
        for jj = 1:k:length(Aa)
            if ii<=jj
                Aa(ii:ii+k-1,jj:jj+k-1) = 0;
            else 
                Aa(ii:ii+k-1,jj:jj+k-1)=A^((1/k)*(ii-jj)-1);
            end
        end
   end
  
   %{
  size(Yy)
  size(c0)
  size(A0)
  size(Xo)
  size(c0*A0*Xo)
  
  size(c0)
  size(Aa)
  size(Bb)
  size(Uu)
  %size(c0*Aa*Bb*Uu)
 %} 
    %Aa(isnan(Aa))=0;
    minimize( norm(Yy-c0*A0*Xo-c0*Aa*Bb*Uu)) %objective function min||cX-Y||
cvx_end
size(Cc);
Cest=Cc;
end
