function evaluation_info=evaluate_KPCSH(XTrain,YTrain,LTrain,XTest,YTest,LTest,param)

    % parameters
    kdim = param.nbits/param.sf;
    
    % normalize
    XTrain = NormalizeFea(XTrain,1);
    YTrain = NormalizeFea(YTrain,1);
    NLTrain = NormalizeFea(LTrain,1);
    XTest = NormalizeFea(XTest,1);
    YTest = NormalizeFea(YTest,1);
    
    keyid = Keypoint(NLTrain, param.thresh);
    % training
    [B,XW,YW] = train_KPCSH(XTrain,YTrain,LTrain,NLTrain,keyid,param);
    
    % evaluate
    BxTest = rsign(XTest*XW,param.nbits,kdim);
    ByTest = rsign(YTest*YW,param.nbits,kdim);
    
    DHamm = BxTest*B';
    [~, orderH] = sort(DHamm,2,'descend');
    evaluation_info.Image_VS_Text_MAP50 = mAP(orderH', LTrain, LTest, 50);
    fprintf('Image_VS_Text_MAP50: %f.\n', evaluation_info.Image_VS_Text_MAP50);
    DHamm = ByTest*B';
    [~, orderH] = sort(DHamm,2,'descend');
    evaluation_info.Text_VS_Image_MAP50 = mAP(orderH', LTrain, LTest, 50);
    fprintf('Text_VS_Image_MAP50: %f.\n', evaluation_info.Text_VS_Image_MAP50);
    
end