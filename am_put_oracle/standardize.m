function y = standardize(x)
n = size(x,1);
d = size(x,2);

if ( d == 1 )
    y = (x - repmat(25, n, d)) ./ repmat(40 - 25, n, d);
end

if ( d == 2 )
    y = (x - repmat(25, n, d)) ./ repmat(30, n, d);
end
    
end