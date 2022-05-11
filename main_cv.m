function main_cv(cv_flag,DATA,prop)
intMat = DATA.interaction;
N = length(prop);  %表示测试样本的比例
option = [3.0000    0.100  100.0000   10.0000];
cv = 5;
result = [];

seed = [1:1000:20000];
if cv_flag==1
    for i = 1:N
        i
        for j = 1:length(seed)
            cv_data = cross_validation(intMat,cv_flag,cv,prop(i),seed(j));
            result0 = five_cross(DATA,cv_data,option);
            result0
            result = [result;[option,prop(i),result0]]; %参数，prop，结果
        end
    end
    best_WHAIMC_CVa1 = result;
    save best_WHAIMC_CVa1 best_WHAIMC_CVa1;
elseif cv_flag==2
    for j = 1:length(seed)
        cv_data = cross_validation(intMat,cv_flag,cv,1,seed(j));
        result0 = five_cross(DATA,cv_data,option);
        result0
        result = [result;[option,result0]]; %参数，结果
    end
    best_WHAIMC_CVr1 = result;
    save best_WHAIMC_CVr1 best_WHAIMC_CVr1;
elseif cv_flag==3
    for j = 1:length(seed)
        cv_data = cross_validation(intMat,cv_flag,cv,1,seed(j));
        result0 = five_cross(DATA,cv_data,option);
        result0
        result = [result;[option,result0]]; %参数，结果
    end
    best_WHAIMC_CVc1 = result;
    save best_WHAIMC_CVc1 best_WHAIMC_CVc1;
end


end

function result = five_cross(DATA,cv_data,option)
intMat = DATA.interaction;
cv = length(cv_data);
result = 0;
for k=1:cv
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%
    W = cv_data{k,1}; %Mask for training set
    train_set = intMat.*W;    %train matrix
%     Y = train_set;
    Y = Assist_method_Utils.WKNKN(train_set,DATA.vir_sim,DATA.drug_sim{1},5,0.9);
    %%%%%%计算病毒相似矩阵(并利用高斯核融合）
    KSNS_Svir = Assist_method_Utils.KSNS_opt(Y);     KSNS_Svir(1:size(KSNS_Svir,1)+1:end)=1;
    K1(:,:,1) =  DATA.vir_sim;  K1(:,:,2) =  KSNS_Svir;
    K1 = Assist_method_Utils.complete_sim(K1);
    [weight_v1] = Assist_method_Utils.cka_kernels_weights(K1,train_set,1);
    vir_sim = Assist_method_Utils.combine_kernels(weight_v1, K1);   
    
    %%%计算药物的相似矩阵
    KSNS_Sdrug = Assist_method_Utils.KSNS_opt(Y');    KSNS_Sdrug(1:size(KSNS_Sdrug,1)+1:end)=1;
    K2(:,:,1) =  DATA.drug_sim{1};K2(:,:,2) =  DATA.drug_sim{2};
    K2(:,:,3) =  KSNS_Sdrug;
    K2 = Assist_method_Utils.complete_sim(K2);
    [weight_v2] = Assist_method_Utils.cka_kernels_weights(K2,train_set,2);
    drug_sim = Assist_method_Utils.combine_kernels(weight_v2,K2);
    
%     W = ones(size(W));
    scores = WHAIMC_opt(W,train_set, vir_sim,drug_sim,option);
    %Result Evaluation    
    pscores = scores(cv_data{k,2});
    result0 = evaluate_opt1(pscores,cv_data{k,3});
    result = result + result0;   
    result/k
    toc;      
end
 result = result/cv;
end

%intMat: interaction matrix
%cv=5,five fold cross
%prop: Ratio of positive and negative cases in the test set
function cv_data = cross_validation(intMat,cv_flag,cv,prop,seed)
cv_data = cell(cv,3);  
[num_A,num_B] = size(intMat);
[pos_indx,pos_indy] = find(intMat==1);
num_pos = length(pos_indx);
if exist('seed','var') == 1
    rand('state',seed)
end
if cv_flag == 1    %interaction cross
    index = randperm(num_pos)';
elseif cv_flag == 2  %row cross
    index = randperm(num_A)';
elseif cv_flag == 3  %col cross
    index = randperm(num_B)';
end
num_step = floor(length(index)/cv);
for j = 1:cv   %
    if j < cv
        ii = index((j-1)*num_step+1:j*num_step);
    else
        ii = index((j-1)*num_step+1:end);
    end
    if cv_flag == 1    %interaction cross
        pos_test_index = [pos_indx(ii),pos_indy(ii)];   %Positive example of test set index
        [neg_indx,neg_indy] = find(intMat==0); %
        num_neg = length(neg_indx);
        %The number of negative samples is prop times of positive samples
        num_test_pos = size(pos_test_index,1);  %The number of positive test samples
        num_test_neg = prop*num_test_pos; %The number of negative test samples

        index00 = randperm(num_neg);
        test_index = [pos_test_index;[neg_indx(index00(1:num_test_neg)),neg_indy(index00(1:num_test_neg))]];
        Test_matrix = full(sparse(test_index(:,1),test_index(:,2),...
            ones(length(test_index(:,1)),1),num_A,num_B));
        test_index = find(Test_matrix==1);    %test single index
        test_label = intMat(test_index);  %test label
        W = ones(num_A,num_B);
        W(test_index) = 0;  %W: train mask matrix
        cv_data(j,:) = {W, test_index, test_label};
    elseif cv_flag == 2   %for all row
        test_index=[];
        for k=1:length(ii)
            test_index = [test_index;[ii(k)*ones(num_B,1),[1:num_B]']];
        end
        %%only one choice
        Test_matrix = full(sparse(test_index(:,1),test_index(:,2),...
            ones(length(test_index(:,1)),1),num_A,num_B));
        test_index = find(Test_matrix==1);
        test_label = intMat(test_index);
        W = ones(num_A,num_B);
        W(test_index) = 0;
        gca = pcolor([W(find(sum(W,2)==0),:);W(find(sum(W,2)==num_B),:)]);
        set(gca, 'LineStyle','none');
        colorbar;
        cv_data(j,:) = {W, test_index, test_label};
        
    elseif cv_flag==3     %对所有列交叉
        test_index=[];
        for k = 1:length(ii)
            test_index = [test_index;[[1:num_A]',ii(k)*ones(num_A,1)]];
        end
        %%only one choice
        Test_matrix = full(sparse(test_index(:,1),test_index(:,2),ones(length(test_index(:,1)),1),num_A,num_B));
        test_index = find(Test_matrix==1);
        test_label = intMat(test_index);
        W = ones(num_A,num_B);
        W(test_index) = 0;
        gca = pcolor([W(:,find(sum(W)==0)),W(:,find(sum(W)==num_A))]);
        set(gca, 'LineStyle','none');
        colorbar;
        cv_data(j,:) = {W, test_index, test_label};
    end
end
end

function result = evaluate_opt1(score,label)
%%%%计算TN,TP,FN,FP
sort_predict_score=unique(sort(score));
score_num = length(sort_predict_score);
Nstep = 2000;
threshold=sort_predict_score(ceil(score_num*(1:Nstep)/(Nstep+1)));

threshold=threshold';
threshold = threshold(end:-1:1);
threshold_num=length(threshold);
TN=zeros(threshold_num,1);
TP=zeros(threshold_num,1);
FN=zeros(threshold_num,1);
FP=zeros(threshold_num,1);

for i=1:threshold_num
    tp_index=logical(score>=threshold(i) & label==1);
    TP(i,1)=sum(tp_index);
    
    tn_index=logical(score<threshold(i) & label==0);
    TN(i,1)=sum(tn_index);
    
    fp_index=logical(score>=threshold(i) & label==0);
    FP(i,1)=sum(fp_index);
    
    fn_index=logical(score<threshold(i) & label==1);
    FN(i,1)=sum(fn_index);
end


%%%%%计算AUPR
SEN=TP./(TP+FN);
PRECISION=TP./(TP+FP);
recall=SEN;
x=recall;
y=PRECISION;
[x,index]=sort(x);
y=y(index,:);
x = [0;x];  y = [1;y];
x(end+1,1)=1;  y(end+1,1)=0;
AUPR=0.5*x(1)*(1+y(1));
for i=1:threshold_num
    AUPR=AUPR+(y(i)+y(i+1))*(x(i+1)-x(i))/2;
end
AUPR_xy = [x,y];

%%%%%计算AUC
AUC_x = FP./(TN+FP);      %FPR
AUC_y = TP./(TP+FN);      %tpR
[AUC_x,ind] = sort(AUC_x);
AUC_y = AUC_y(ind);
AUC_x = [0;AUC_x];
AUC_y = [0;AUC_y];
AUC_x = [AUC_x;1];
AUC_y = [AUC_y;1];

AUC = 0.5*AUC_x(1)*AUC_y(1)+sum((AUC_x(2:end)-AUC_x(1:end-1)).*(AUC_y(2:end)+AUC_y(1:end-1))/2);

AUCxy = [AUC_x(:),AUC_y(:)];

%%%%%计算其他指标
temp_accuracy=(TN+TP)/length(label);   %%准确率
temp_sen=TP./(TP+FN);    %%真实的有链接，预测正确的正确率
recall = temp_sen;
temp_spec=TN./(TN+FP);   %%真实无连接，预测正确率
temp_precision=TP./(TP+FP); %%预测有链接的正确率
temp_f1=2*recall.*temp_precision./(recall+temp_precision);
[~,index]=max(temp_f1);
%%%%计算F1最高的如下值：
f1=temp_f1(index);
%%%%计算得分前10，前15，前20的召回率
[~,ind] = sort(score,'descend');

% precision_top5 = sum(label(ind(1:5)))/sum(label);
% precision_top10 = sum(label(ind(1:10)))/sum(label);
% precision_top20 = sum(label(ind(1:20)))/sum(label);
% precision_top30 = sum(label(ind(1:30)))/sum(label);
% precision_top40 = sum(label(ind(1:40)))/sum(label);
prop = round(length(label)*[0.02,0.04,0.06,0.08,0.1]);
precision_top_p2 = sum(label(ind(1:prop(1))))/sum(label);
precision_top_p4 = sum(label(ind(1:prop(2))))/sum(label);
precision_top_p6 = sum(label(ind(1:prop(3))))/sum(label);
precision_top_p8 = sum(label(ind(1:prop(4))))/sum(label);
precision_top_p10 = sum(label(ind(1:prop(5))))/sum(label);


result=[AUPR,AUC,f1,...
    precision_top_p2,precision_top_p4, precision_top_p6,precision_top_p8,...
    precision_top_p10];
%     precision_top10,precision_top20, precision_top30,precision_top40,...

end
