function [M,sca_M,fp,syn,rsqv,stf_out,G0]=moment_invt(green,obs,lenstf,lambda,sublen,option)
%==========================================================================
% moment_invt is a function for seismic moment inversion. it can retrieve the
% seismic moment tensor and source time function at the same time
%--------------------------------------------------------------------------
% Input:
%      green: the green's functions for all stations. it has a zise of [leng,6,nsta].
%             where leng is the length of data for sampling dots. 6 is the number of
%             elements of moment tensor. nsta is the number of stations.
%      obs:   size of [leng,nsta]
%     lenstf: length of source time function for sampling dots
%     lambda: weight of minimum seismic moment tensor constraint
%     sublen: ==0: smooth the source time function
%             ==1: do not smooth the source time function
%             ==2,3,4...:length of sub-stf which consist the source time function
%     option: ==1: do not constrain seismic moment tensor in any way (full moment) 
%             ==2: only clear isotropic component
%          others: clear isotropic component and CLVD component    
% Output:
%        M: solution for seismic moment tensor
%     sca_M: scalar seismic moment
%       fp: fault parameters for the best double couple
%      rsqv: residual for each interation
%     stf_out: the source time function
%
%--------------------------------------------------------------------------
%                             Zhang Yong, Xu Lisheng
%                             2008/09/04
%==========================================================================
%==========================================================================
%����ע�ͣ�
%moment_invt�������������ݵ���������ĺ��������ܹ�ͬʱ�õ�����غ���Դʱ�亯����
%--------------------------------------------------------------------------
%�� �룺
%  green������̨վ�ĸ��ֺ��������ľ����С��[leng,6,nsta]�����leng�����ǲ���
%         �����ݵĳ��ȡ�6�Ǿ����������ĸ�����nsta��̨վ�ĸ�����
%  obs��  ����[leng,nsta]�Ĵ�С
% lenstf����������Դʱ�亯���ĳ���
% lambda����С�ĵ���ؾ�����Ȩ��Լ��
% sublen��==0��������Դʱ�亯�����⻬
%         ==1��  �����⻬
%         ==2,3,4...������Դʱ�亯�������ĳ��ȣ����ǹ������ܵ�ʱ�亯��
% option��==1��  ���Ե�����������κ�Լ������������أ�
%         ==2��  ֻȥ������ͬ�Գɷ�
%      others��  ȥ������ͬ�Ժ�CLVD�ɷ�
%
%�� ����
%       M������������Ľ�
%   sca_M�����������
%      fp�����˫��ż�Ķϲ����
%    rsqv��ÿ�ε����Ĳв�
% stf_out����Դʱ�亯��
%
%--------------------------------------------------------------------------
%                              ���£�������
%                              2008/09/04
%==========================================================================
%
%
% the default inversion is for best double couple source
% Ĭ�Ϸ��ݵ������˫��żԴ
if nargin<6; option=3;
    % default sets for source time function is smoothing
    % Ĭ�ϵ������Ƕ���Դʱ�亯�����⻬
    if nargin<5; sublen=0;
        % if no neccessary, we do not constrain the scalar moment
        % ���û�б�Ҫ�Ļ���һ�㲻�Ա�������Լ��
        if nargin<4; lambda=0;
        end
    end
end

% size of green's functions
% ���ֺ����ľ����С
sg=size(green);

% maximum interation times. usually 5 interations is enough
% ������������һ��5�ε������㹻��
itern=8;

% rep==0 while mod(lenstf,sublen)==0
% ��lengstf/sublen��������0����rep==0
rep=0;
if sublen>0&&mod(lenstf,sublen)~=0
    % while mod(lenstf,sublen)~=0, we set rep to be 1
    % ��lengstf/sublen����������0������rep��ֵ��1
    rep=1;
    % preserve the origin length of source time function 
    % ������Դʱ�亯����ԭʼ����
    lenstf0=lenstf;
    % correct the length of source time function insure that mod(lenstf,sublen)==0 
    % У����Դʱ�亯���ĳ��ȣ��Ա�֤mod(lenstf,sublen)==0
    lenstf=ceil(lenstf./sublen)*sublen;
end

% initialization of source time function for the interation
% ���ڵ�������Դʱ�亯���ĳ�ʼ��������ֵ
stf_in=ones(lenstf./max(sublen,1),1);

% residual of interations
% �����Ĳв�
rsqv=zeros(itern,1);
% elements of seismic moment tensors in all interations
% ���е����еĵ���صĸ�����
M=zeros(6,itern);
% source time functions in all interations 
% ���дε�������Դʱ�亯��
stf=zeros(length(stf_in),itern);


%-------------------------------------------------------------------------
if sublen<0
    error('length of triangles must be positive number!');
end
if sublen>0
    
    % triangles consisting source time function
    % ����ʱ�亯����������
    substf=tri_stf(sublen,sublen);

    % green0 is the convolution of green and substf 
    % green0�Ǹ��ֺ�������ʱ�亯���ľ��
    green0=green;
    for i=1:sg(2)
        for j=1:sg(3)
            zg=conv(green(:,i,j),substf(2:end));
            green0(:,i,j)=zg(1:sg(1));
        end
    end
else
    green0=green;
end
%-------------------------------------------------------------------------

% the interation 
% ����
for kk=1:itern

    if sublen>1
        % convert stf_in0 with a sampling rate of 1/sublen to stf_in1 with a sampling rate of 1
        % ��������Ϊ1/sublen��stf_in0ת��Ϊ������Ϊ1��stf_in1
        stf_in0=interm(stf_in,[sublen,1]);
%        stf_in1=[[0:stf_in0(1)/(sublen-1):stf_in0(1)]';stf_in0(2:end)];
        stf_in1=[[0:stf_in0(1)/(sublen):stf_in0(1)]';stf_in0(2:end)];
    else
        stf_in1=stf_in;
    end

    % convolution of green and source time function
    % ���ֺ�������Դʱ�亯���ľ��
    g=zeros(sg(1)+length(stf_in1)-1,sg(2),sg(3));
    for i=1:sg(3)
        for j=1:sg(2)
            g(:,j,i)=conv(green(:,j,i),stf_in1);
        end
    end
    g(sg(1)+1:end,:,:)=[];

    % construction of green's matrix
    % ���ֺ�������Ĺ���
    G=zeros(sg(1)*sg(3),sg(2));
    for i=1:sg(3)
        G((i-1)*sg(1)+1:i*sg(1),:)=g(:,:,i);
    end

    % solve the seismic moment tensor 
    % �����������
    if option==1
        % for full seismic moment tensor invertion
        % �����ĵ������������
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500);
    elseif option==2
        % only set the trace of moment tensor to be 0
        % ֻ��Ҫ���õ������������ļ�Ϊ0
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500,'minX=sum(X([1,4,6]))./3;X([1,4,6])=X([1,4,6])-minX;');
    else
        % only preserve DC (Double Couple) components
        % ֻ����DC��˫��ż���ɷ�
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500,'minX=sum(X([1,4,6]))./3;X([1,4,6])=X([1,4,6])-minX;[m,fp]=mt2fp0_zh(X);X=fp2mt_zh(m,fp(1,1),fp(1,2),fp(1,3));');
    end
    
    % get the residual in this interation
    % �õ��ôε����Ĳв�
    rsqv(kk)=min(rsq);

    g0=zeros(sg(1),sg(3));
    for i=1:sg(3)
        for j=1:6
            g0(:,i)=g0(:,i)+green0(:,j,i)*M(j,kk);
        end
    end

    % matrix for source time function inversion
    % ��Դʱ�亯�����ݵľ���
    G0=zeros(sg(1)*sg(3),length(stf_in));
    for i=1:sg(3)
        subG=con_subG_M(g0(:,i),lenstf,max(sublen,1),sg(1));
        G0((i-1)*sg(1)+1:i*sg(1),:)=subG;
    end

    sG=size(G0);
    lam=mean(abs(G0(:)));
    
    if sublen==0
        % smooth the source time function 
        % ��Դʱ�亯���Ĺ⻬
        %stf_in=cgls_xu0([G0;ones(1,sG(2))*lam*lambda],[obs(:);0],500,'X([1,end])=0;X(X<0)=0;X=sm(X,3);X([1,end])=0;');
        stf_in=cgls_xu0([G0;ones(1,sG(2))*lam*lambda],[obs(:);0],500,'X([1,end])=0;X(X<0)=0;X=sm(X,3);X([1,end])=0;');
    elseif sublen==1
        % sublen==1 
        stf_in=cgls_xu0([G0;ones(1,sG(2))*lam*lambda],[obs(:);0],500,'X([1,end])=0;X(X<0)=0;');
    else
        % sublen>1 
        stf_in=cgls_xu0([G0;ones(1,sG(2))*lam*lambda],[obs(:);0],500,'X([end])=0;X(X<0)=0;');
    end

    % preserve the source time function
    % ������Դʱ�亯��
    stf(:,kk)=stf_in;
end

% find the result with minimum residual
% Ѱ�Ҿ�����С�в�Ľ�
[m,n]=min(rsqv);
M=M(:,n);

if n==1
    % for interation number equals 1
    % ���ڵ�����������1�����
    stf_out=stf_in;
else
    % for interation number more than 1
    % ���ڵ�����������1�����
    stf_out=stf(:,n-1);
end

if sublen>1
    % resample the source time function
    % ����Դʱ�亯���ز���
    stf_out=interm(stf_out,[sublen,1]);
    stf_out=[[0:stf_out(1)/(sublen-1):stf_out(1)]';stf_out(2:end)];
end

% calculate the fault parameters from
% seismic moment tensor
% �ӵ������������ϲ����
[sca_M,fp]=mt2fp0_zh(M);
sca_M=sca_M*sum(stf_out);
syn=reshape(G*M,size(obs));

% while mod(lenstf,sublen)~=0, clear the
% dots complemented before
% ��mod(lenstf,sublen)~=0ʱ�����֮ǰ�����ĵ�
if rep==1;
    stf_out(lenstf0+1:end)=[];
end
return
%=================================end======================================


function subG=con_subG_M(gstf,lensub,d_time,leng)
%

if mod(lensub,d_time)~=0
    error('lensub and d_time are not fit');
end

subG=zeros(leng,lensub./d_time);

len=length(gstf);
for i=1:lensub./d_time
    %    subG((i-1)*d_time+[1:len],i)=gstf;
    subG((i-1)*d_time+1:min((i-1)*d_time+len,end),i)=gstf(1:min(leng-(i-1)*d_time,end));
    %    subG((i)*d_time+1:min((i)*d_time+len,end),i)=gstf(1:min(leng-(i)*d_time,end));
end

function stf = tri_stf(len1,len2)
% =================================================================================
% stf = tri_stf(len1,len2)
%     TRI_STF is to generate a triang function with rising time of len1 and
%     descending time of len2�� Get the result as stf ,with the length
%     len1+len2+1
% -------------------------------------------------------------------------
%                                    Zhang Yong Xu Lisheng and Chen Yuntai
%                                             2006/10/14/15:50  ����������
% ==================================================================================
% ==================================================================================
% stf = tri_stf(len1,len2)
%     TRI_STF����������һ�������κ�������������ʱ����len1���½�ʱ����len2���õ�
%     �����stf��������len1+len2+1
% -------------------------------------------------------------------------
%                                      ���£�������������̩
%                                    2006/10/14/15:50   ����������
% ==================================================================================
stf=[0:1/len1:1,1-1/len2:-1/len2:0];
stf=stf(:);
return

function b=interm(a,x)
%b=interm(a,x)
sa=size(a);
if isscalar(x)
    b0=interp1(0:sa(1)-1,a,0:1/x:sa(1)-1);
    b00=interp1(0:sa(2)-1,b0',0:1/x:sa(2)-1);
    b=b00';
else
    b0=interp1(0:sa(1)-1,a,0:1/x(1):sa(1)-1);
    if sa(2)~=1
        b00=interp1(0:sa(2)-1,b0',0:1/x(2):sa(2)-1);
        b0=b00;
    end
    b=b0';
end

function [X,rsq]=cgls_xu0(A,b,itern,enchar)
%---------initial
%��ʼ��
sa=size(A);
X=zeros(sa(2),1);
X0=zeros(sa(2),itern);
s0=b-A*X; %s0: xi
g=A'*s0; %r0,p0
h=g; %r0,p0
rsq=zeros(itern,1);
%-------------------------------
for iter=1:itern
    q=A*h; %xi=A*h
    anum=sum(g.*h);
    aden=q'*q;
    if aden==0
        disp('very singular matrix')
    end
    anum=anum./aden;

    X=X+anum.*h;
    %constrains: positive constraints------------------------------
    %Լ������ֵԼ��
    %     xabs=abs(X);
    if nargin>3
        eval(enchar);
    end

    X0(:,iter)=X;
    %--------------------------------------------------------------
    s=(A*X)-b;%xj=(A*X)-b;
    %s=s0-anum*q;
    rsq(iter)=s'*s;
    %rp=rsq;
    r=A'*s;%xi=A'*xj;
    gg=sum(g.^2);
    dgg=sum((r+g).*r);
    g=-r;
    h=g+dgg./gg.*h;
end
[m,n]=min(rsq);
X=X0(:,n);
return
%--------------------------End----------
%------------------------------

function a=sm(b,n)
% si=size(b);
% if si(1)<si(2)
%     b=b';
% end
si=size(b);
a=zeros(si);
if n==1
    a=b;
    return
end
if mod(n,2)==0
    nx=n/2;
    %bsm=[zeros(nx,si(2));b;zeros(nx,si(2))];
    bsm=[ones(nx,si(2))*diag(b(1,:));b;ones(nx,si(2))*diag(b(end,:))];
    for i=1:si(1)
        a(i,:)=mean(bsm(i:i+n-1,:));
    end
    %     for j=1:si(2)
    %         bsm=[zeros(nx,1);b(:,j);zeros(nx-1,1)];
    %         for i=1:si(1)
    %             a(i,j)=mean(bsm(i:i+n-1));
    %         end
    %     end
else
    nx=(n-1)/2;
    %bsm=[zeros(nx,si(2));b;zeros(nx,si(2))];
    bsm=[ones(nx,si(2))*diag(b(1,:));b;ones(nx,si(2))*diag(b(end,:))];
    for i=1:si(1)
        a(i,:)=mean(bsm(i:i+n-1,:));
    end
    %     for j=1:si(2)
    %         bsm=[zeros(nx,1);b(:,j);zeros(nx,1)];
    %         for i=1:si(1)
    %             a(i,j)=mean(bsm(i:i+n-1));
    %         end
    %     end
end

function [sca_M,fp]=mt2fp0_zh(M)
% =========================================================================
%  [sca_M,fp]=mt2fp0_zh(M);
%  MT2FP0_ZH is a function to solve focal mechanisms for moment tensors
% -------------------------------------------------------------------------
% Input:
%        M: a matrix with the size of [6,num]. where num is the number of
%           moment tensor.
%           for each seismic moment tensor, we only use the 6 independent
%           elements. for the i-th moment tensor, as an example:
%           M(:,i)=[M(1,1),M(1,2),M(1,3),M(2,2),M(2,3),M(3,3)]'
%
% Output:
%   sca_M: the scalar moment. size: [num,1]
%      fp: the mechanisms. size: [num*2,3]
%          every two rows correspond the two plane for a mechanism.
%          the first column is strike, secondary is dip, last is rake
%
%  Notice: this function is developed based on mt2fp of Xu, it can run fast
%          while solve large number of moment tensors at the same time.
%          here the traces of each moment tensor should be set to be 0s, or
%          else there may be errors or warnings while calculating.
%
% for more details,please see mt2fp
% -------------------------------------------------------------------------
%                                 Zhang Yong, Chen Yun-Tai, Xu Li-Sheng
%                                  2007/01/26/01:37 ,  Peking University
% =========================================================================
% =========================================================================
%   [sac_M,fp]=mt2fp0_zh(M);
%  MT2FP0_ZH������������ʽ����Դ���ƽ�ĺ���
% -------------------------------------------------------------------------
% �� �룺
%      M����СΪ[6,num]�ľ�������num�Ǿ��������ܸ�����
%         ����ÿһ�����������������ֻ��Ҫʹ��6�������ķ��������ڵ�i���������
%         ����������˵��
%         M(:,i)=[M(1,1),M(1,2),M(1,3),M(2,2),M(2,3),M(3,3)]'
%
% �� ����
%  sca_M����������أ���С��[num,1]
%     fp����Դ���ƣ���С��[num*2,3]
%         ÿ���ж�Ӧ��һ����Դ���Ƶ��������档��һ�������򣬵ڶ�������ǣ�������
%         �ǻ����ǡ�
%
%  ע �⣺�˺����ǻ���Xu��mt2fp�Ŀ�������ͬʱ����������������ʱ����������
%         ��ĸ��졣���ÿһ������������ļ���Ӧ������Ϊ0����������ڼ����ʱ
%         ����־�������
%
% �����Ҫ�˽����ϸ�ڣ���鿴mt2fp
% -------------------------------------------------------------------------
%                                        ���£�����̩��������
%                                    2007/01/26/01:37��������ѧ
% =========================================================================
sM=size(M);
nx=[1,2,3,2,4,5,3,5,6];

if nargout==1
    sca_M=zeros(sM(2),1);
    for i=1:sM(2)
        mom=reshape(M(nx,i),[3,3]);
        [eigval]=eig(mom);
        %eigval=diag(eigval);

        sca_M(i)=(max(eigval(:))-min(eigval(:)))./2;
    end
    return
else
    eigvec0=zeros(3,3*sM(2));
    eigval0=zeros(3,3*sM(2));
    %fp=zeros(2*sM(2),3);
    addvec=3:3:3*sM(2)+2;
    for i=1:sM(2)
        mom=reshape(M(nx,i),[3,3]);
        numth=addvec(i);
        [eigvec0(:,numth-2:numth),eigval0(:,numth-2:numth)]=eig(mom);
        %eigval=diag(eigval0(numth-2:numth,:));
    end
    eigval0=reshape(eigval0,[9,sM(2)]);
    % only keep 1,5,9--the diag elements, delete others:
    % ֻ��Ҫ����1,5,9�������Խ���Ԫ�أ�Ȼ��ɾ��������
    eigval0([2,3,4,6,7,8],:)=[];
    [mt,nt]=max(eigval0);
    [mp,np]=min(eigval0);
    sca_M=(mt-mp)./2;
    sca_M=sca_M';

    % find positions of tvec and pvec
    % �ҵ�tvec��pvec��λ��
    nt=nt+addvec-3;
    np=np+addvec-3;% here addvec-3=0:3:3*sM(2)-1;

    % select all tvecs and pvecs, and get two matrixs
    % ѡ�����е�tvec��pvec��Ȼ��õ���������
    tvec=eigvec0(:,nt);
    pvec=eigvec0(:,np);

    % Plane I:
    % ����I�� 
    uone=tvec+pvec;
    vone=tvec-pvec;
    % Plane II
    % ����II��
    utwo=vone;
    vtwo=uone;

    %check Plane I:
    %if vone(3)>0,vone=-vone;uone=-uone;end
    %��֤����I�����vone(3)>0,����vone=-vone;uone=-uone;Ȼ�����
    v13g0=vone(3,:)>0;
    vone(:,v13g0)=-vone(:,v13g0);
    uone(:,v13g0)=-uone(:,v13g0);
    %check Plane II:
    %if vtwo(3)>0,vtwo=-vtwo;utwo=-utwo;end
    %��֤����II�����vtwo(3)>0,����vtwo=-vtwo;utwo=-utwo;Ȼ�����
    v23g0=vtwo(3,:)>0;
    vtwo(:,v23g0)=-vtwo(:,v23g0);
    utwo(:,v23g0)=-utwo(:,v23g0);

    %set the u-v matrix to be vector, and put together
    %����u-v�����Ϊ������Ȼ�����Ƿ���һ��
    u=sqrt(2)./2.*reshape([uone;utwo],[3,sM(2)*2]);
    v=sqrt(2)./2.*reshape([vone;vtwo],[3,sM(2)*2]);

    dip=acos(-v(3,:));
    strike=atan2(-v(1,:),v(2,:));
    rake=atan2(-u(3,:)./sin(dip),u(1,:).*cos(strike)+u(2,:).*sin(strike));

    %Strikes are agreed to be positive
    %����ͨ������Ϊ����ֵ
    strikel0=find(strike<0);
    strike(strikel0)=strike(strikel0)+2*pi;
    fp=[strike' dip' rake'];

    %convert to Degree
    %������תΪ�Ƕ�
    fp=fp./pi.*180;
    return
end
% ==================================end=
% ===================================


function moment=fp2mt_zh(m,x,y,z)
%moment=momentt(m,x,y,z)
%m: scalar seismic moment
%x: strike
%y: dip
%z: rake
%moment: seismic moment
m=m(:)';
x=x(:)';
y=y(:)';
z=z(:)';
moment=zeros(6,length(m));

x=x.*pi./180;
y=y.*pi./180;
z=z.*pi./180;

sinx=sin(x);cosx=cos(x);sin2x=2.*sinx.*cosx;cos2x=cosx.*cosx-sinx.*sinx;
siny=sin(y);cosy=cos(y);sin2y=2.*siny.*cosy;cos2y=cosy.*cosy-siny.*siny;
sinz=sin(z);cosz=cos(z);
moment(1,:)=-m.*(cosz.*siny.*sin2x+sinz.*sin2y.*(sinx.*sinx));
moment(2,:)=m.*(cosz.*siny.*cos2x+sinz.*sin2y.*sin2x./2);
%moment(2,1)=moment(1,2);
moment(3,:)=-m.*(cosz.*cosy.*cosx+sinz.*cos2y.*sinx);
%moment(3,1)=moment(1,3);
moment(4,:)=m.*(cosz.*siny.*sin2x-sinz.*sin2y.*cosx.*cosx);
moment(5,:)=m.*(-cosz.*cosy.*sinx+sinz.*cos2y.*cosx);
%moment(3,2)=moment(2,3);
moment(6,:)=m.*sinz.*sin2y;
return
% ================================end======================================
% moment(1,:)=-m.*(cos(z).*sin(y).*sin(2.*x)+sin(z).*sin(2.*y).*(sin(x)).^2);
% moment(2,:)=m.*(cos(z).*sin(y).*cos(2.*x)+(1./2).*sin(z).*sin(2.*y).*sin(2.*x));
% %moment(2,1)=moment(1,2);
% moment(3,:)=-m.*(cos(z).*cos(y).*cos(x)+sin(z).*cos(2.*y).*sin(x));
% %moment(3,1)=moment(1,3);
% moment(4,:)=m.*(cos(z).*sin(y).*sin(2.*x)-sin(z).*sin(2.*y).*(cos(x)).^2);
% moment(5,:)=m.*(-cos(z).*cos(y).*sin(x)+sin(z).*cos(2.*y).*cos(x));
% %moment(3,2)=moment(2,3);
% moment(6,:)=m.*sin(z).*sin(2.*y);