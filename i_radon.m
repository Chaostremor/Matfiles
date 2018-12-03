function b = i_radon(m, dt, x, q, flag)
% inv_radon ��Radon�򷵻ص�t-x�� ����ΪRadon���е����ݾ���ÿ��q��Ϊһ�С�
%
% ���� :
%     m : Radon���е����� ( m[nt, nq], nt��ʱ�����������nq��q������Ŀ��
%     dt��ʱ������������ sΪ��λ
%     x��ƫ�ƾ�������λ�ã���mΪ��λ 
%     q�����߲��� ��flag=1�������������� �� flag=2��
%     flag=1: tau-p �任�� flag=2:  tau-q�任
% �����
%     b : ���ص�x-t ���еĵ������� ��b[nt, nx]��

[nt, nq] =  size(m);   % ��Radon�����ݵ�ά��
nx = length(x);   % �������ĵ���
fp = fft(m, [], 1);   % tau-p/q��任��f-q/p��

R = zeros(nq, nx);
temp = zeros(nx, 1);
ufxx = zeros(nt, nx);

for k = 1:nt   % f-q/p��任��f-x��
    if k <= nt/2+1
        omega=6.28318530717959*(k-1)/(nt*dt);
        R=exp(i*omega*q'*(x.^flag));   % �������Ի������� L����
        % R=exp(i*omega*q'*((x.^2+4*q.^2).^0.5-2*q));   % ����˫���� L����
        temp = R'*reshape(fp(k,1:nq),nq,1); 
        ufxx(k,:) = temp.';
    else
        ufxx(k,:) = conj(ufxx(nt+2-k,:));
    end
end

b=real(ifft(ufxx,[],1));%f-x��任��t-x��
