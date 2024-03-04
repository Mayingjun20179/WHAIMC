classdef Assist_method_Utils
    methods (Static)
        %%%%%%%%%%   KSNS
        function S = KSNS_opt(X,sim)
            if nargin==1
                distance = pdist(X);
                distance = squareform(distance);  %计算欧式距离
                sim = max(distance(:))-distance;
            end
            neighbor_num = floor(0.3*size(sim,1));
            feature_matrix = X;
            nearst_neighbor_matrix = Assist_method_Utils.KSNS_neighbors(sim,neighbor_num);
            S = Assist_method_Utils.jisuanW_opt(feature_matrix,nearst_neighbor_matrix);
        end
        
        function nearst_neighbor_matrix=KSNS_neighbors(sim,neighbor_num)
            %%%%nearst_neighbor_matrix：  represent Neighbor matrix
            N = size(sim,1);
            D = sim-diag(inf*gpuArray.ones(1,N));
            [~, si]=sort(D,2,'descend');
            nearst_neighbor_matrix = gpuArray.zeros(N,N);
            index=si(:,1:neighbor_num);
            for i=1:N
                nearst_neighbor_matrix(i,index(i,:))=1;
            end
        end
        
        %The iterative process of this algorithm
        function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix)
            lata1 = 4;  lata2 = 1;
            X = feature_matrix';  % each column of X represents a sample, and each behavior is an indicator
            [~,N] = size(X);    % N represents the number of samples
            C = nearst_neighbor_matrix';
            rand('state',1);
            W = gpuArray(single(rand(N,N)));
            W = W- diag(diag(W));
            W = W./repmat(sum(W),N,1);
            G  = Assist_method_Utils.jisuan_Kel(X);
            G(isnan(G))=0;
            G = G/max(G(:));
            
            WC1 = W'*G*W-2*W*G+G;
            WC = sum(diag(WC1))/2;
            wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
            objv = wucha;
            jingdu = 0.0001;
            error = jingdu*(1+lata1+lata2);   %Iteration error threshold
            we = 1;      %Initial value of error
            gen=1;
            while  gen<100 && we>error
                %gen
                FZ = G+lata1*C.*W;
                FM = G*W+lata1*W+lata2*W;
                FM(FM==0) = eps;  W = FZ./FM.*W;
                %     W = FZ./(FM+eps).*W;
                WC1 = W'*G*W-2*W*G+G;
                WC = sum(diag(WC1))/2;
                objv1 = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
                we = abs(objv1-objv(end));
                objv = [objv,objv1];
                gen = gen+1;
            end
            W=W';
            W = Assist_method_Utils.matrix_normalize(W);
        end
        
        function W = matrix_normalize(W)
            K = 10;
            W(isnan(W))=0;
            W(1:(size(W,1)+1):end)=0;
            for round=1:K
                SW = sum(W,2);
                ind = find(SW>0);
                SW(ind) = 1./sqrt(SW(ind));
                D1 = diag(SW);
                W=D1*W*D1;
            end
        end
        
        function K  =jisuan_Kel(X)
            %X Columns represent samples, and rows represent features
            X = X';
            sA = (sum(X.^2, 2));
            sB = (sum(X.^2, 2));
            K = exp(bsxfun(@minus,bsxfun(@minus,2*X*X', sA), sB')/mean(sA));
        end
        
        
        
        %%%空缺相似性矩阵整合
        function K = complete_sim(K)
            W = size(K,3);
            for i=1:W
                [M,N] = size(K(:,:,i));
                for j = 1:M
                    for k=1:N
                        if K(j,k,i) == 0
                            K(j,k,i) = max(K(j,k));
                        end
                    end
                end
            end
        end
        %%%%多网络整合
        function [w] = cka_kernels_weights(Kernels_list,adjmat,dim)
            
            % adjmat : binary adjacency matrix
            % dim    : dimension (1 - rows, 2 - cols)
            
            num_kernels = size(Kernels_list,3);
            weight_v = zeros(num_kernels,1);
            y = adjmat;
            % Graph based kernel
            if dim == 1
                ga = y*y';
            else
                ga = y'*y;
            end
            %ga = Knormalized(ga);
            N_U = size(y,dim);
            l=ones(N_U,1);
            U = eye(N_U) - (l*l')/N_U;   %中心矩阵
            M = zeros(num_kernels,num_kernels);
            for i=1:num_kernels
                for j=1:num_kernels
                    kk1 = U*Kernels_list(:,:,i)*U;
                    kk2 = U*Kernels_list(:,:,j)*U;
                    mm = trace(kk1'*kk2);
                    m1 = trace(kk1*kk1');
                    m2 = trace(kk2*kk2');
                    M(i,j) = mm/(sqrt(m1)*sqrt(m2));
                end
            end
            a = zeros(num_kernels,1);
            for i=1:num_kernels
                kk = U*Kernels_list(:,:,i)*U;
                aa = trace(kk'*ga);
                a(i) = aa*(N_U-1)^-2;
            end
            v = randn(num_kernels,1);
            falpha = @(v)Assist_method_Utils.obj_function(v,M,a);
            [x_alpha, fval_alpha] = Assist_method_Utils.optimize_weights(v, falpha);
            w = x_alpha;
        end
        
        function [J] = obj_function(v,M,a)
            J =v'*M*v - v'*a ;
        end
        
        function [x, fval] = optimize_weights(x0, fun)
            n = length(x0);
            Aineq   = [];
            bineq   = [];
            Aeq     = ones(1,n);
            beq     = 1;
            LB      = zeros(1,n);
            UB      = ones(1,n);
            options = optimoptions('fmincon','Algorithm','interior-point', 'Display', 'notify');
            [x,fval] = fmincon(fun,x0,Aineq,bineq,Aeq,beq,LB,UB,[],options);
        end
        %根据权重对多核进行整合
        function result = combine_kernels(weights, kernels)
            % length of weights should be equal to length of matrices
            n = length(weights);
            result = zeros(size(kernels(:,:,1)));
            for i=1:n
                result = result + weights(i) * kernels(:,:,i);
            end
        end
        
        
        
        function Y_tem = WKNKN(Y, SD, ST, K, eta)
            Yd = zeros(size(Y));
            Yt = Yd;
            wi = zeros(K,1);
            wj = zeros(1,K);
            [num_drugs, num_targets] = size(Y);
            for i = 1:num_drugs
                [~,ind] = sort(SD(i,:),'descend');
                dnn_i = ind(2:K+1);
                Zd = sum(SD(i, dnn_i));
                for ii = 1:K
                    wi(ii) = eta^ii* SD(i,dnn_i(ii));
                end
                if abs(Zd)>eps
                    Yd(i,:) = sum(repmat(wi,1,num_targets).*Y(dnn_i,:)) / Zd;
                end
            end
            
            for j = 1:num_targets
                [~,ind] = sort(ST(j,:),'descend');
                tnn_j = ind(2:K+1);
                Zt = sum(ST(j, tnn_j));
                for jj = 1:K
                    wj(jj) = eta^jj * ST(j,tnn_j(jj));
                end
                if abs(Zt)>eps
                    Yt(:,j) = sum(repmat(wj,num_drugs,1).*Y(:,tnn_j),2) / Zt;
                end
            end
            
            Ydt = (Yd + Yt)/2;
            Y_tem = max(Ydt,Y);
        end
        
    end
end