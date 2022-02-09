function v = rpdf(X, Y, x, y, r, s)
    v = exp(-0.5*((sqrt((X-x).^2+(Y-y).^2)-r)./s).^2);
%     if sum(v) < numel(v)*1e-10
%         v = ones(size(v));
%     end
end