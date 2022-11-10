function Zmne = mne(G2dU_red, Ca, lambda)
        [Nch, Nsites_red] = size(G2dU_red);
        Nsites_red = Nsites_red / 2;
        Cs = eye(Nsites_red * 2);
        Cn = eye(size(G2dU_red, 1));
        GGt = G2dU_red * Cs * G2dU_red';
        Wmne = G2dU_red' * inv(G2dU_red * Cs * G2dU_red' + lambda * trace(GGt) / size(GGt, 1) * Cn);

        clear Zmne;
        range2d = 1:2;
        for i = 1:Nsites_red
            w = Wmne(range2d, :);
            m = w * Ca * w';
            [u ss v] = svd(m);
            Zmne(i) = ss(1, 1);
            range2d = range2d + 2;
        end
end
