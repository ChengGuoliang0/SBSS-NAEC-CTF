function [e, W] = SBSS_NAEC_CTF_AuxIVA(x, y, p, L, nfft, alpha)

Lx=length(x);
xp=x.^(2*(1:p)-1);
yxp=[y,xp];

win=hanning(nfft,'periodic');
shift=fix(nfft/4);
N=fix((Lx+nfft)/shift)-1;
nf=fix(nfft/2)+1;
X=zeros(p+1,N,nf);
for i=1:p+1
    X(i,:,:)=StFT(yxp(:,i),nfft,win,shift).';
end

beta=0.4;
Xn=zeros(p*L+1,1,nf);
W=repmat([1,zeros(1,p*L)],1,1,nf);
Sn=zeros(nf,1);
V=repmat(1e-3*eye(p*L+1),1,1,nf);
e1=[1;zeros(p*L,1)];
E=zeros(nf,N);

fprintf('%d:     ',N);
for n=1:N
    fprintf('\b\b\b\b\b%5d', n);
    Xn=[X(:,n,:);Xn(2:p*(L-1)+1,:,:)];
    for k=1:nf
        Sn(k)=W(:,:,k)*Xn(:,:,k);
    end
    r=sqrt(sum(abs(Sn).^2));
    phi=(r+eps).^(beta-2);
    for k=1:nf
        V(:,:,k)=alpha*V(:,:,k)+(1-alpha)*phi*Xn(:,:,k)*Xn(:,:,k)';
        W(:,:,k)=((V(:,:,k)+1e-6*eye(p*L+1))\e1)';
        W(:,:,k)=W(:,:,k)/W(1,1,k);
        E(k,n)=W(:,:,k)*Xn(:,:,k);
    end
end
e=iStFT(E,Lx,win,shift);

end