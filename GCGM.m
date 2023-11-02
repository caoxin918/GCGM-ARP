function[x_k1] = GCGM(Nodes,b_analytic,G,alpha)

NumNodes = size(Nodes, 1);
%初始化变量
x_k = pinv(G)*b_analytic;
z_k = zeros(NumNodes, 1);
alpha = 1e-5;
lambda = 10;
%最大迭代次数
k_max = 100;
% beta = 0.02;
beta = 1e-3;
A = G;
A_inv = G';
y = b_analytic;
% E_L2_sum = 0;
% E_cos_sum = 0;

% N = 20;
for k = 1:k_max
    % 确定下降方向 z_k
    temp1 = beta / lambda
    t1 = (1 - temp1) * x_k;
    t2 = A_inv * (A * x_k - y) /lambda;
%     t1 = (1 + beta / (lambda * norm(x_k))) * x_k;
%     t2 = A_inv * (A * x_k - y) /lambda;
    t = t1 - t2;
    % t = t / max(abs(t));
    temp = alpha / lambda;
    
    % idx1 = find(t >= temp);
    % idx2 = find(abs(t) < temp);
    % idx3 = find(t <= -temp);
    % z_k(idx1,1) = t(idx1,1) - temp;
    % z_k(idx2,1) = 0;
    % z_k(idx3,1) = t(idx3,1) + temp;
    idx1 = find(t >= temp);
    idx2 = find(t < -temp);
    idx3 = find(t > -temp & t  < temp);
    z_k(idx1,1) = sign(t(idx1)) .* (abs(t(idx1)) - temp);
    z_k(idx2,1) = sign(t(idx1)) .* (abs(t(idx1)) + temp);
    z_k(idx3,1) = 0;
%     确定步长 sk
    s_k = determinS(x_k,z_k);
%     更新x_k1
    x_k1 = x_k + s_k * (z_k - x_k);
    x_k = x_k1;
    
    % E_L2 = norm(G * x_k - b_analytic);
    % E_cos = sum(G * x_k .* b_analytic) / (norm(G * x_k) * norm(b_analytic));
    % 
    % E_L2_sum = E_L2_sum + E_L2;
    % E_cos_sum = E_cos_sum + E_cos;
    % 
    % P_L2 = (1 ./ E_L2_sum) / sum(1 ./ E_L2_sum);
    % P_cos = E_cos_sum / sum(E_cos_sum);
    % P_err = (P_L2 + P_cos) / 2;
    % P_X = x_k .* P_err / sum(x_k .* P_err);
   
end



