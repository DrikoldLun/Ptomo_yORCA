function id_stanw = stanw_unique(sta,nwk)
    for i = 1:length(sta)
        stanwk{i} = [sta{i},nwk{i}];
    end
    stanwk = stanwk';
    stanwkuni = unique(stanwk);
    id_stanw = zeros(length(stanwk),1);
    for i = 1:length(stanwkuni)
        id = find(strcmp(stanwk,stanwkuni{i}));
        id_stanw(i) = id(1);
    end
end