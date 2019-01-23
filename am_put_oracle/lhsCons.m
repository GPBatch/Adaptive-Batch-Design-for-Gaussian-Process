% Generate lhs design in a constrained space
function y = lhsCons(len, rec)

d = size(rec,1);
y = lhsdesign(len,d);

start = repmat(rec(:,1)',len,1);
ed = repmat(rec(:,2)',len,1);
diff = ed - start;
y = y.*diff + start;
end

