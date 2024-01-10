function M = contact_list(N,u,dir)
   C_Lis = rangesearch(N(:,1:3),N(:,1:3),(3)*max(N(:,4)));
    M = zeros(4*size(N,1),1);
    ll=1;

    for ii=1:1:size(N,1)
        B = N(C_Lis{ii,1},:);
        for jj=2:1:size(B,1)
            Nxx = N(ii,1) - B(jj,1) ;
            Nyy = N(ii,2) - B(jj,2) ;
            Nzz = N(ii,3) - B(jj,3) ;
            r_eff = 2*N(ii,4)*B(jj,4)/(N(ii,4)+B(jj,4));
            D = sqrt(Nxx^2+Nyy^2+Nzz^2);
            Ndd = D - N(ii,4) - B(jj,4);
            if (Ndd <= u*r_eff && N(ii,8)< B(jj,8))%%%Considering the effect of h
                M(ll,1:4) = [N(ii,8),B(jj,8),D,asind(abs(N(ii,dir)-B(jj,dir))/D)];
                ll = ll+1;
            end 
        end
    end
end

