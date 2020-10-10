function Zp = Apsiicos(G2dU_red, Ca, Upwr, Rnk, pos_correction, type, Wps)

% pos_correction: 0 for changing for zeros and 1 for absolute values

    if strcmp(type,'power')
        P = Upwr(:,1:Rnk)*Upwr(:,1:Rnk)';
    elseif strcmp(type,'whitened')
        UcorrRnk = Upwr(:,1:Rnk);
        P = inv(Wps)*(eye(size(UcorrRnk,1))-UcorrRnk*UcorrRnk')*Wps;
    else
        P = eye(Nch^2)-Upwr(:,1:Rnk)*Upwr(:,1:Rnk)';
    end
    
    [Nch, Nsites_red] = size(G2dU_red);
    Nsites_red = Nsites_red/2;
    Cap = reshape(P*Ca(:),size(Ca));            
    
    [e a] = eig(Cap);
    if pos_correction == 1;
        Cap = e*abs(a)*e';
        iCap = tihinv(Cap, 0.001);
    else 
        a(a<0) = 0;
        Cap = e*abs(a)*e';
        iCap = tihinv(Cap, 0.05);
    end
    
    range2d = 1:2;
    for i = 1:Nsites_red
        g = G2dU_red(:,range2d);
        m = inv(g'*iCap*g);
        [u ss v] = svd(m);
        Zp(i) = ss(1,1);
        range2d = range2d+2;
    end
end