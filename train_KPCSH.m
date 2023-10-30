function [B,XW,YW] = train_KPCSH(XTrain,YTrain,LTrain,NLTrain,keyid,param)
    
    % parameters
    max_iter = param.iter;
    kdim = param.nbits/param.sf;
    omega = param.omega;
    theta = param.theta;
    beta = param.beta;
    alpha = param.alpha;
    lambda = param.lambda;
    mu = param.mu;
    nbits = param.nbits;
    sel_num = 1000;
    if strcmp(param.db_name, 'NUS-WIDE')
        sel_num = 10000;
    end
    
    n = size(LTrain,1);
    m = size(keyid,1);
    dX = size(XTrain,2);
    dY = size(YTrain,2);
    nonkeyid = setdiff(1:n,keyid);
    NLTrain1 = NLTrain(keyid,:);
    NLTrain2 = NLTrain(nonkeyid,:);

    % hash code learning
    if kdim < n
        H = sqrt(m*nbits/kdim)*orth(rand(m,kdim));
        B1 = rsign(H,nbits,kdim);
        B2 = rsign(sqrt((n-m)*nbits/kdim)*orth(rand((n-m),kdim)),nbits,kdim);
        for i = 1:max_iter
            % update H
            Z = omega * B1 + nbits* NLTrain1*(NLTrain1'*B1) + nbits *alpha* NLTrain1 * (NLTrain2' * B2);
            [~,Lmd,VV] = svd(Z'*Z);
            index = (diag(Lmd)>1e-6);
            V = VV(:,index); V_ = orth(VV(:,~index));
            U = Z *  (V / (sqrt(Lmd(index,index))));
            U_ = orth(randn(m,kdim-length(find(index==1))));
            H = sqrt(m*nbits/kdim)*[U U_]*[V V_]';

            % update B1
            B1 = rsign(omega * H + nbits * NLTrain1*(NLTrain1'*H)+theta*m*nbits/kdim*ones(m,kdim),nbits,kdim);
            
            % update B2
            B2 = rsign(nbits * alpha * NLTrain2*(NLTrain1'*H)+theta*(n-m)*nbits/kdim*ones(n-m,kdim),nbits,kdim);
        end
    end
    clear Z Temp Lmd VV index U U_ V V_ BCluster
    
    B = zeros(n,kdim);
    B(keyid,:) = B1;
    B(nonkeyid, :) = B2;
    
    save KPCSHcode B keyid
    

    % hash function learning
    if m < sel_num
        sel_idx = keyid;
    else
        sel_idx = randperm(m,sel_num);
    end
    Bs = B(sel_idx,:);
    YW = rand(dY,kdim);
    
    if mu == 0
        max_iter = 1;
    end
    for i = 1:max_iter
        XW = (XTrain'*XTrain+(lambda)*eye(dX))\(XTrain'*B+ mu*XTrain'*YTrain*YW +((XTrain'*NLTrain)*NLTrain(sel_idx,:)')*Bs*beta*nbits)...
        /((1+mu) * eye(kdim)+Bs'*Bs*beta);
        YW = (YTrain'*YTrain+lambda*eye(dY))\(YTrain'*B+ mu*YTrain'*XTrain*XW +((YTrain'*NLTrain)*NLTrain(sel_idx,:)')*Bs*beta*nbits)...
        /((1+mu) * eye(kdim)+Bs'*Bs*beta);
    end
end