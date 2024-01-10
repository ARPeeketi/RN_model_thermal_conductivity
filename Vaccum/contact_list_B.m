function M = contact_list_B(A,B,u)
        ll=1;
        M=[];
     for ii =1:1:size(A,1)
         for jj = 1:1:size(B,1)
            Nxx = A(ii,1) - B(jj,1) ;
            Nyy = A(ii,2) - B(jj,2) ;
            Nzz = A(ii,3) - B(jj,3) ;
            r1 = A(ii,4);
            r2 = B(jj,4);
            r_eff = 2*r1*r2/(r1+r2);
            Ndd = sqrt(Nxx^2+Nyy^2+Nzz^2) - r1 - r2;
            if Ndd < u*r_eff  %%%Considering the effect of h
                M(ll,1) = A(ii,8);
                M(ll,2) = B(jj,8);
                M(ll,3) = sqrt(Nxx^2+Nyy^2+Nzz^2);
                M(ll,4) = asind(abs(Nzz)/sqrt(Nxx^2+Nyy^2+Nzz^2));
                ll = ll+1;
            end 
         end
    end

end