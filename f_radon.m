% could refer to paper <高精度Radon变换方法及去噪分析，郭全仕>，
% <Radon变换的MALAB实现，沈操>
% Author: , modified: C. Song
function m = f_radon(b , dt , x , q , flag ,noise)
% 正Radon变换。输入地震数据，返回Radon域中的数据
% b: 地震数据（b【nt ， nx】， nt是时间采样点数，nx是道数，即每列数据是一道）
% dt：时间采样间隔， 以s为单位，
% x：偏移距或地震道的位置，以m为单位
% q：射线参数（flag=1）或者抛物线曲率（flag=2）
% flag：1 线性，tau-p变换；2 抛物线， tau-q变换.
% noise：白噪系数
% m：radon域中的数据
[nt,nx]=size(b);
nq=length(q);
uf=fft(b,[],1);
R=zeros(nq,nx);
toep=zeros(nq,nx);
g=zeros(nq,1);
ufpp=zeros(nt,nq);
%f-x域变换到f-p域
for k=1:nt %对每一个频率成分   线性抛物线
    %由于是函数的fourier变换是hermit函数，故采用共轭对称赋值
    if k <= nt/2+1
        omega=6.28318530717959*(k-1)/(nt*dt);
        R=exp(i*omega*q'*(x.^flag));%计算线性或抛物线 L矩阵
        % R=exp(i*omega*q'*((x.^2+4*q.^2).^0.5-2*q));%计算双曲线 L矩阵
        
        MATRIX=R*R';
        toep=inv(MATRIX+noise*eye(nq))*R;      % 阻尼最小二乘方法求广义逆
        g=toep*reshape(uf(k,1:nx),nx,1);
        
        ufpp(k,:)=g.';
    else
        ufpp(k,:)=conj(ufpp(nt+2-k,:));
    end
end

m=real(ifft(ufpp,[],1));%f-p域变换到tau-q/p域
