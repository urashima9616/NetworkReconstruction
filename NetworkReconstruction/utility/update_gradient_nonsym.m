function new_gradient = update_gradient_nonsym(N, A_obs, P_est, P_k_est, labels, cur_gradient, idx_j, idx_k, psi, psi_new)
new_gradient = cur_gradient;
for i=1:N
    %Update the contribution of edge (i, idx_j)
    if A_obs(i, idx_j) == 1
        %Remove the contribution of this edge
        mapped_i = psi(i);
        mapped_j = psi(idx_j);
        new_mapped_j = psi_new(idx_j);
        
        i_labels = labels(mapped_i,:);
        j_labels = labels(mapped_j,:);
        new_j_labels = labels(new_mapped_j, :);
        
        temp_cnt = zeros(size(cur_gradient));
        temp_cnt_new = temp_cnt;
        %Get the profile of label hittings
        for kk = 1:length(i_labels)
            temp_cnt(i_labels(kk), j_labels(kk)) =  temp_cnt(i_labels(kk), j_labels(kk)) + 1;
            temp_cnt_new(i_labels(kk), new_j_labels(kk)) =  temp_cnt_new(i_labels(kk), new_j_labels(kk)) + 1;
        end
        
        for l = 1:length(cur_gradient(1,:))
            for n = 1 : length(cur_gradient(1,:))
                %(P_uv/1-P_uv + 1)*(t/P_i,j)
                old = (temp_cnt(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, mapped_j))^(-1));
                new = (temp_cnt_new(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, new_mapped_j))^(-1));
                diff = new - old;
                
                new_gradient(l,n) = new_gradient(l,n) + diff;
                %new_gradient(n,l) = new_gradient(l,n);
                %gradient(l,n) = gradient(l,n) - 2*(temp_cnt(l,n)/P_est(l,n));
            end
        end
    end
    
    
    
    %Update the contribution of edge (i, idx_k)
    if A_obs(i, idx_k) == 1 && i ~= idx_j % edge (idx_k, idx_j) considered above already
        %Remove the contribution of this edge
        mapped_i = psi(i);
        mapped_k = psi(idx_k);
        new_mapped_k = psi_new(idx_k);
        
        i_labels = labels(mapped_i,:);
        k_labels = labels(mapped_k,:);
        new_k_labels = labels(new_mapped_k, :);
        
        temp_cnt = zeros(size(cur_gradient));
        temp_cnt_new = temp_cnt;
        %Get the profile of label hittings
        for kk = 1:length(i_labels)
            temp_cnt(i_labels(kk), k_labels(kk)) =  temp_cnt(i_labels(kk), k_labels(kk)) + 1;
            temp_cnt_new(i_labels(kk), new_k_labels(kk)) =  temp_cnt_new(i_labels(kk), new_k_labels(kk)) + 1;
        end
        
        for l = 1:length(cur_gradient(1,:))
            for n = 1 : length(cur_gradient(1,:))
                %(P_uv/1-P_uv + 1)*(t/P_i,j)
                old = (temp_cnt(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, mapped_k))^(-1));
                new = (temp_cnt_new(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, new_mapped_k))^(-1));
                diff = new - old;
                
                new_gradient(l,n) = new_gradient(l,n) + diff;
                %new_gradient(n,l) = new_gradient(l,n);
                %gradient(l,n) = gradient(l,n) - 2*(temp_cnt(l,n)/P_est(l,n));
            end
        end
    end
    
end
