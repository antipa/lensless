function A = generate_2D_MURA(N)
%N =3329;
ismod = @(q,n)any(~mod((0:n).^2-q,n));
An = zeros(N,1);
for n = 0:N-1

       ncheck = ismod(n,N);
        if ~ncheck
            cn = -1;
        else
            cn = 1;
        end
       
        An(n+1) = cn;

end
A = double(An*An'>0);
A(:,1) = 1;
A(1,:) = 0;
imagesc(A)