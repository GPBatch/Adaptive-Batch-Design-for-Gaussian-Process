function EI = cf_metric(cand_mean, cand_sd, nug, design)

switch design
    case 'MCU'
        % calculate the local empirical error %
        EI = normcdf(-abs(cand_mean)./cand_sd);
    case 'cSUR'
        sigmanoise = nug;
        % approximate the new one-step-ahead variance %
        cand_varnew = sigmanoise .* cand_sd.^2./(sigmanoise + cand_sd.^2);
        % calculate the reduced empirical error %
        EI = normcdf( -abs(cand_mean)./cand_sd ) + normcdf( -abs(cand_mean)./sqrt(abs(cand_varnew)) );
end
end