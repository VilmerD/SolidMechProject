function dwdzf = sCAeigs(B, phi, w, K, M, Kold, R, phi0, edof, np, K0, ...
    M0, feig)
% Generic method, i don't ork explain now

if strcmp(feig, 'cae')
    dwdzf = sCAE(B, phi, w, K, M, Kold, R, phi0, edof, np, K0, M0);
elseif strcmp(feig, 'caeeon')
    dwdzf = sCAEEON(B, phi, w, K, M, Kold, R, phi0, edof, np, K0, M0);
elseif strcmp(feig, 'caeeonmod')
    dwdzf = sCAEEONmod(B, phi, w, K, M, Kold, R, phi0, edof, np, K0, M0);
elseif strcmp(feig, 'caers')
    dwdzf = sCAERS(B, phi, w, K, M, Kold, R, phi0, edof, np, K0, M0);
elseif strcmp(feig, 'caeeonlazy')
    dwdzf = Deigen(phi, w, edof, K0, M0);
else
    % It should not be possible to select wrong basis generation method and
    % end up here, since this method is always called after CAeigs which
    % has the same check.
    msg = sprintf('Adjoint method matching "%s" not found', feig);
    errorStruct.message = msg;
    error(errorStruct);
end

end