% could refer to paper <�߾���Radon�任������ȥ���������ȫ��>��
% <Radon�任��MALABʵ�֣����>
% Author: , modified: C. Song
function m = f_radon(b , dt , x , q , flag ,noise)
% ��Radon�任������������ݣ�����Radon���е�����
% b: �������ݣ�b��nt �� nx���� nt��ʱ�����������nx�ǵ�������ÿ��������һ����
% dt��ʱ���������� ��sΪ��λ��
% x��ƫ�ƾ��������λ�ã���mΪ��λ
% q�����߲�����flag=1���������������ʣ�flag=2��
% flag��1 ���ԣ�tau-p�任��2 �����ߣ� tau-q�任.
% noise������ϵ��
% m��radon���е�����
[nt,nx]=size(b);
nq=length(q);
uf=fft(b,[],1);
R=zeros(nq,nx);
toep=zeros(nq,nx);
g=zeros(nq,1);
ufpp=zeros(nt,nq);
%f-x��任��f-p��
for k=1:nt %��ÿһ��Ƶ�ʳɷ�   ����������
    %�����Ǻ�����fourier�任��hermit�������ʲ��ù���ԳƸ�ֵ
    if k <= nt/2+1
        omega=6.28318530717959*(k-1)/(nt*dt);
        R=exp(i*omega*q'*(x.^flag));%�������Ի������� L����
        % R=exp(i*omega*q'*((x.^2+4*q.^2).^0.5-2*q));%����˫���� L����
        
        MATRIX=R*R';
        toep=inv(MATRIX+noise*eye(nq))*R;      % ������С���˷����������
        g=toep*reshape(uf(k,1:nx),nx,1);
        
        ufpp(k,:)=g.';
    else
        ufpp(k,:)=conj(ufpp(nt+2-k,:));
    end
end

m=real(ifft(ufpp,[],1));%f-p��任��tau-q/p��
