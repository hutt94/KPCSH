function keypointlist = Keypoint(L_tr, thresh)
    [a,~,c]=unique(L_tr,'row');
    for i = 1:size(a,1)
        count(i)=sum(c==i);
    end
    a = NormalizeFea(a,1);
    S = a*a';
    degree = sum(S,1)./sum(S>0,1);
    key = (count/max(count)).*degree;

    [~,b]=sort(key,'descend');

    keyid = b(1:thresh);

    keypointlist = [];
    for i = 1: thresh
        id = find(c==keyid(i));
        keypointlist = [keypointlist;id];
    end
end