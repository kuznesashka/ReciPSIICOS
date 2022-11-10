function [Upwr, ds, A] = ProjectorOnlyAwayFromPowerComplete(G2dU, PwrRnk, fixed)

if(nargin < 2)
    PwrRnk = 350;
end

if fixed == 0
    Ns = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % perform projection of the coherence matrix away from the power only
    Nch = size(G2dU, 1);

    % project for each potential source
    fprintf('Collecting power subspace of coherence span...\n');
    A = zeros(Nch^2, Ns * 3);
    for i = 1:Ns
         gi = G2dU(:, 2 * i - 1);
         v = gi * gi';
         A(:, 3 * i - 2) = v(:) / norm(v(:));
         gj = G2dU(:, 2 * i);
         v = gj * gj';
         A(:, 3 * i - 1) = v(:) / norm(v(:));
         v = gi * gj' + gj * gi';
         A(:, 3 * i) = v(:) / norm(v(:));
    end
else
    Ns = size(G2dU, 2);
    Nch = size(G2dU, 1);

    % project for each potential source
    fprintf('Collecting power subspace of coherence span...\n');
    A = zeros(Nch ^ 2, Ns);
    for i = 1:Ns
         gi = G2dU(:, i);
         v = gi * gi';
         A(:, i) = v(:) / norm(v(:));
    end
end

fprintf('Finding eigen space...\n');
AA = A * A';
[u s] = eigs(AA, PwrRnk);
ds = diag(s);
Upwr = u(:, 1:PwrRnk);

end
