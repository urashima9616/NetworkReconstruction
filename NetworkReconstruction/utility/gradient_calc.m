%This is a utility to calculate the gradient the LL(theta)
%
%const1 = -k||P_est||_1^(k-1)
%const2 = -k(||P_est||_2^(2))^(k-1)*P_est_{i,j}

function gradient = gradient_calc(P_est, P_k_est, m, k,  G_sample, psi)

gradient = zeros(size(P_est));

%Create map between edge and a counting matrix
%keySet = {'00'};
%valueSet = {zeros(size(P_est))};
%mapObj = containers.Map(keySet,valueSet);

%Get the edgelist from the graph
e = G_sample.Edges;
edgelist = table2array(e);

%Iterate on the edge list
for i = 1:length(edgelist(:,1))
    %Get the edge name
    %edge_name = sprintf('$d%d', edgelist(i,1), edgelist(i,2));
    %tf = isKey(mapObj, edge_name);
    
    %If it exists
    mapped_i = psi( edgelist(i,1));
    mapped_j = psi( edgelist(i,2));
    
    
    i_labels = calculate_labels(m, k, mapped_i);
    j_labels = calculate_labels(m, k, mapped_j);
    
%     if tf == 0
%         mapObj(edge_name) = zeros(size(P_est));
%     end
    
    temp_cnt = zeros(size(gradient));
    
    %Get the profile of label hittings 
    for kk = 1:length(i_labels)
        temp_cnt(i_labels(kk), j_labels(kk)) =  temp_cnt(i_labels(kk), j_labels(kk)) + 1;
        
    end
    %mapObj(edge_name) = temp_cnt;
    
    %Update the gradient matrix with the contribution of current edge
    for l = 1:length(gradient(1,:))
        for n = l : length(gradient(1,:))
            %(P_uv/1-P_uv + 1)*(t/P_i,j)
            gradient(l,n) = gradient(l,n) + (temp_cnt(l,n)/P_est(l,n))*((1-P_k_est(mapped_i, mapped_j))^(-1));
            %gradient(l,n) = gradient(l,n) - 2*(temp_cnt(l,n)/P_est(l,n));
        end
    end
    
  
end
    
    

%Calculate the # of P_est_{i,j} in each edge

const1 = -k*sum(sum(P_est))^(k-1);
const2 = -k*(sum(sum(P_est.^(2))))^(k-1);
const3 = -k*(sum(sum(P_est.^(3))))^(k-1);

for l = 1:length(gradient(1,:))
    for n = l:length(gradient(1,:))
         gradient(l,n) = gradient(l,n) + const1 + const2*P_est(l,n) + const3*(P_est(l,n)^(2));   
         gradient(n,l) = gradient(l,n);
    end
end
