function s = sample_pmf(v, w, n)
% Draw from pmf
    edges = cumsum([0; w(:)]);
    edges = edges/edges(end);
%     if isnan(sum(edges))
%         hai
%     end
    [~, ~, i] = histcounts(rand(1, n), edges);
    s = v(i, :);
end