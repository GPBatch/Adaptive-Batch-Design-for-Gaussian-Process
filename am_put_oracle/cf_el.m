%%%%%%%%%%%%%%%
%' Expected Loss for Contour Finding
%'
%'  Compute expected loss using the optimal stopping loss function.
%' @param obj  must be a DynaTree or a list with an sd field
%' @export
function el = cf_el(objMean,objSd)

el = objSd.*normpdf( -abs(objMean)./objSd ) - abs(objMean).*normcdf( -abs(objMean)./objSd);
el(el<0) = 0;

end