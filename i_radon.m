function b = i_radon(m, dt, x, q, flag)
% inv_radon 从Radon域返回到t-x域。 输入为Radon域中的数据矩阵，每个q道为一列。
%
% 输入 :
%     m : Radon域中的数据 ( m[nt, nq], nt是时间采样点数，nq是q道的数目）
%     dt：时间采样间隔，以 s为单位
%     x：偏移距地震道的位置，以m为单位 
%     q：射线参数 （flag=1）或抛物线曲率 （ flag=2）
%     flag=1: tau-p 变换； flag=2:  tau-q变换
% 输出：
%     b : 返回到x-t 域中的地震数据 （b[nt, nx]）

[nt, nq] =  size(m);   % 求Radon域数据的维数
nx = length(x);   % 求地震道的道数
fp = fft(m, [], 1);   % tau-p/q域变换到f-q/p域

R = zeros(nq, nx);
temp = zeros(nx, 1);
ufxx = zeros(nt, nx);

for k = 1:nt   % f-q/p域变换到f-x域
    if k <= nt/2+1
        omega=6.28318530717959*(k-1)/(nt*dt);
        R=exp(i*omega*q'*(x.^flag));   % 计算线性或抛物线 L矩阵
        % R=exp(i*omega*q'*((x.^2+4*q.^2).^0.5-2*q));   % 计算双曲线 L矩阵
        temp = R'*reshape(fp(k,1:nq),nq,1); 
        ufxx(k,:) = temp.';
    else
        ufxx(k,:) = conj(ufxx(nt+2-k,:));
    end
end

b=real(ifft(ufxx,[],1));%f-x域变换到t-x域
