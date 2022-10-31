function [e, W] = SBSS_NAEC_CTF_ILRMA(x, y, p, L, nb, nfft, alpha)

Lx=length(x);
xp=x.^(2*(1:p)-1);
yxp=[y,xp];
win=hanning(nfft,'periodic');
shift=fix(nfft/4);
J=fix((Lx+nfft)/shift)-1;
I=fix(nfft/2)+1;
X=zeros(p+1,J,I);
for i=1:p+1
    X(i,:,:)=StFT(yxp(:,i),nfft,win,shift).';
end
N=p*L+1;
W=repmat([1,zeros(1,N-1)],1,1,I);
T = max( rand( I, nb ), eps );
Y = zeros(I,1);
Xj = zeros(N,1,I);
D = repmat(eye(N),1,1,I);
e1=[1;zeros(N-1,1)];
E = zeros(I,J);

fprintf('%d:     ',J);
for j=1:J
    fprintf('\b\b\b\b\b%5d', j);
    Xj=[X(:,j,:);Xj(2:p*(L-1)+1,:,:)];
    V = max( rand( nb, 1 ), eps );
    R = T*V;
    for i=1:I
        Y(i) = W(:,:,i)*Xj(:,:,i);
    end
    P = max(abs(Y).^2,eps);
    T = T .* sqrt( (P.*(R.^(-2)))*V.' ./ ( (R.^(-1))*V.' ) );
    T = max(T,eps);
    R = T*V;
    V = V .* sqrt( T.'*(P.*(R.^(-2))) ./ ( T.'*(R.^(-1)) ) );
    V = max(V,eps);
    R = T*V;
    for i=1:I
        D(:,:,i) = alpha*D(:,:,i)+(1-alpha)*((Xj(:,:,i)./R(i,:))*Xj(:,:,i)');
        W(:,:,i) = ((D(:,:,i)+1e-6*eye(N))\e1)';
        W(:,:,i)=W(:,:,i)/W(1,1,i);
        E(i,j) = W(:,:,i)*Xj(:,:,i);
    end
end
e=iStFT(E,Lx,win,shift);

end