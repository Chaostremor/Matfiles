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
%中文注释：
%moment_invt代码是用来反演地震矩张量的函数。它能够同时得到地震矩和震源时间函数。
%--------------------------------------------------------------------------
%输 入：
%  green：所有台站的格林函数。它的矩阵大小是[leng,6,nsta]。这里，leng变量是采样
%         点数据的长度。6是矩张量分量的个数。nsta是台站的个数。
%  obs：  矩阵[leng,nsta]的大小
% lenstf：采样点震源时间函数的长度
% lambda：最小的地震矩矩张量权重约束
% sublen：==0：　对震源时间函数做光滑
%         ==1：  不做光滑
%         ==2,3,4...：子震源时间函数函数的长度，他们构成了总的时间函数
% option：==1：  不对地震矩张量做任何约束（完整地震矩）
%         ==2：  只去掉各向同性成分
%      others：  去掉各向同性和CLVD成分
%
%输 出：
%       M：地震矩张量的解
%   sca_M：标量地震矩
%      fp：最佳双力偶的断层参数
%    rsqv：每次迭代的残差
% stf_out：震源时间函数
%
%--------------------------------------------------------------------------
%                              张勇，许力生
%                              2008/09/04
%==========================================================================
%
%
% the default inversion is for best double couple source
% 默认反演的是最佳双力偶源
if nargin<6; option=3;
    % default sets for source time function is smoothing
    % 默认的设置是对震源时间函数做光滑
    if nargin<5; sublen=0;
        % if no neccessary, we do not constrain the scalar moment
        % 如果没有必要的话，一般不对标量距做约束
        if nargin<4; lambda=0;
        end
    end
end

% size of green's functions
% 格林函数的矩阵大小
sg=size(green);

% maximum interation times. usually 5 interations is enough
% 最大迭代次数，一般5次迭代就足够了
itern=8;

% rep==0 while mod(lenstf,sublen)==0
% 当lengstf/sublen的余数是0，则rep==0
rep=0;
if sublen>0&&mod(lenstf,sublen)~=0
    % while mod(lenstf,sublen)~=0, we set rep to be 1
    % 当lengstf/sublen的余数不是0，设置rep的值是1
    rep=1;
    % preserve the origin length of source time function 
    % 保留震源时间函数的原始长度
    lenstf0=lenstf;
    % correct the length of source time function insure that mod(lenstf,sublen)==0 
    % 校正震源时间函数的长度，以保证mod(lenstf,sublen)==0
    lenstf=ceil(lenstf./sublen)*sublen;
end

% initialization of source time function for the interation
% 用于迭代的震源时间函数的初始化，赋初值
stf_in=ones(lenstf./max(sublen,1),1);

% residual of interations
% 迭代的残差
rsqv=zeros(itern,1);
% elements of seismic moment tensors in all interations
% 所有迭代中的地震矩的各分量
M=zeros(6,itern);
% source time functions in all interations 
% 所有次迭代的震源时间函数
stf=zeros(length(stf_in),itern);


%-------------------------------------------------------------------------
if sublen<0
    error('length of triangles must be positive number!');
end
if sublen>0
    
    % triangles consisting source time function
    % 构成时间函数的三角形
    substf=tri_stf(sublen,sublen);

    % green0 is the convolution of green and substf 
    % green0是格林函数和子时间函数的卷积
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
% 迭代
for kk=1:itern

    if sublen>1
        % convert stf_in0 with a sampling rate of 1/sublen to stf_in1 with a sampling rate of 1
        % 将采样率为1/sublen的stf_in0转换为采样率为1的stf_in1
        stf_in0=interm(stf_in,[sublen,1]);
%        stf_in1=[[0:stf_in0(1)/(sublen-1):stf_in0(1)]';stf_in0(2:end)];
        stf_in1=[[0:stf_in0(1)/(sublen):stf_in0(1)]';stf_in0(2:end)];
    else
        stf_in1=stf_in;
    end

    % convolution of green and source time function
    % 格林函数和震源时间函数的卷积
    g=zeros(sg(1)+length(stf_in1)-1,sg(2),sg(3));
    for i=1:sg(3)
        for j=1:sg(2)
            g(:,j,i)=conv(green(:,j,i),stf_in1);
        end
    end
    g(sg(1)+1:end,:,:)=[];

    % construction of green's matrix
    % 格林函数矩阵的构建
    G=zeros(sg(1)*sg(3),sg(2));
    for i=1:sg(3)
        G((i-1)*sg(1)+1:i*sg(1),:)=g(:,:,i);
    end

    % solve the seismic moment tensor 
    % 求解地震矩张量
    if option==1
        % for full seismic moment tensor invertion
        % 完整的地震矩张量反演
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500);
    elseif option==2
        % only set the trace of moment tensor to be 0
        % 只需要设置地震距张量矩阵的迹为0
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500,'minX=sum(X([1,4,6]))./3;X([1,4,6])=X([1,4,6])-minX;');
    else
        % only preserve DC (Double Couple) components
        % 只保留DC（双力偶）成分
        [M(:,kk),rsq]=cgls_xu0(G,obs(:),500,'minX=sum(X([1,4,6]))./3;X([1,4,6])=X([1,4,6])-minX;[m,fp]=mt2fp0_zh(X);X=fp2mt_zh(m,fp(1,1),fp(1,2),fp(1,3));');
    end
    
    % get the residual in this interation
    % 得到该次迭代的残差
    rsqv(kk)=min(rsq);

    g0=zeros(sg(1),sg(3));
    for i=1:sg(3)
        for j=1:6
            g0(:,i)=g0(:,i)+green0(:,j,i)*M(j,kk);
        end
    end

    % matrix for source time function inversion
    % 震源时间函数反演的矩阵
    G0=zeros(sg(1)*sg(3),length(stf_in));
    for i=1:sg(3)
        subG=con_subG_M(g0(:,i),lenstf,max(sublen,1),sg(1));
        G0((i-1)*sg(1)+1:i*sg(1),:)=subG;
    end

    sG=size(G0);
    lam=mean(abs(G0(:)));
    
    if sublen==0
        % smooth the source time function 
        % 震源时间函数的光滑
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
    % 存下震源时间函数
    stf(:,kk)=stf_in;
end

% find the result with minimum residual
% 寻找具有最小残差的解
[m,n]=min(rsqv);
M=M(:,n);

if n==1
    % for interation number equals 1
    % 对于迭代次数等于1的情况
    stf_out=stf_in;
else
    % for interation number more than 1
    % 对于迭代次数大于1的情况
    stf_out=stf(:,n-1);
end

if sublen>1
    % resample the source time function
    % 对震源时间函数重采样
    stf_out=interm(stf_out,[sublen,1]);
    stf_out=[[0:stf_out(1)/(sublen-1):stf_out(1)]';stf_out(2:end)];
end

% calculate the fault parameters from
% seismic moment tensor
% 从地震矩张量计算断层参数
[sca_M,fp]=mt2fp0_zh(M);
sca_M=sca_M*sum(stf_out);
syn=reshape(G*M,size(obs));

% while mod(lenstf,sublen)~=0, clear the
% dots complemented before
% 当mod(lenstf,sublen)~=0时，清除之前增补的点
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
%     descending time of len2。 Get the result as stf ,with the length
%     len1+len2+1
% -------------------------------------------------------------------------
%                                    Zhang Yong Xu Lisheng and Chen Yuntai
%                                             2006/10/14/15:50  地球物理所
% ==================================================================================
% ==================================================================================
% stf = tri_stf(len1,len2)
%     TRI_STF是用来产生一个三角形函数，它的上升时间是len1，下降时间是len2。得到
%     结果是stf，长度是len1+len2+1
% -------------------------------------------------------------------------
%                                      张勇，许力生，陈运泰
%                                    2006/10/14/15:50   地球物理所
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
%初始化
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
    %约束：正值约束
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
%  MT2FP0_ZH是求解矩张量形式的震源机制解的函数
% -------------------------------------------------------------------------
% 输 入：
%      M：大小为[6,num]的矩阵，其中num是矩张量的总个数。
%         对于每一个地震矩张量，我们只需要使用6个独立的分量。对于第i个地震矩张
%         量，举例来说：
%         M(:,i)=[M(1,1),M(1,2),M(1,3),M(2,2),M(2,3),M(3,3)]'
%
% 输 出：
%  sca_M：标量地震矩，大小：[num,1]
%     fp：震源机制，大小：[num*2,3]
%         每两行对应于一个震源机制的两个节面。第一列是走向，第二列是倾角，第三列
%         是滑动角。
%
%  注 意：此函数是基于Xu的mt2fp的开发，当同时求解许多地震矩张量的时候，它可以运
%         算的更快。这里，每一个地震矩张量的迹都应该设置为0，否则可能在计算的时
%         候出现警告或错误。
%
% 如果想要了解更多细节，请查看mt2fp
% -------------------------------------------------------------------------
%                                        张勇，陈运泰，许力生
%                                    2007/01/26/01:37，北京大学
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
    % 只需要保留1,5,9——即对角线元素，然后删掉其他的
    eigval0([2,3,4,6,7,8],:)=[];
    [mt,nt]=max(eigval0);
    [mp,np]=min(eigval0);
    sca_M=(mt-mp)./2;
    sca_M=sca_M';

    % find positions of tvec and pvec
    % 找到tvec和pvec的位置
    nt=nt+addvec-3;
    np=np+addvec-3;% here addvec-3=0:3:3*sM(2)-1;

    % select all tvecs and pvecs, and get two matrixs
    % 选择所有的tvec和pvec，然后得到两个矩阵
    tvec=eigvec0(:,nt);
    pvec=eigvec0(:,np);

    % Plane I:
    % 节面I： 
    uone=tvec+pvec;
    vone=tvec-pvec;
    % Plane II
    % 节面II：
    utwo=vone;
    vtwo=uone;

    %check Plane I:
    %if vone(3)>0,vone=-vone;uone=-uone;end
    %验证节面I：如果vone(3)>0,则做vone=-vone;uone=-uone;然后结束
    v13g0=vone(3,:)>0;
    vone(:,v13g0)=-vone(:,v13g0);
    uone(:,v13g0)=-uone(:,v13g0);
    %check Plane II:
    %if vtwo(3)>0,vtwo=-vtwo;utwo=-utwo;end
    %验证节面II：如果vtwo(3)>0,则做vtwo=-vtwo;utwo=-utwo;然后结束
    v23g0=vtwo(3,:)>0;
    vtwo(:,v23g0)=-vtwo(:,v23g0);
    utwo(:,v23g0)=-utwo(:,v23g0);

    %set the u-v matrix to be vector, and put together
    %设置u-v矩阵成为向量，然后将它们放在一起
    u=sqrt(2)./2.*reshape([uone;utwo],[3,sM(2)*2]);
    v=sqrt(2)./2.*reshape([vone;vtwo],[3,sM(2)*2]);

    dip=acos(-v(3,:));
    strike=atan2(-v(1,:),v(2,:));
    rake=atan2(-u(3,:)./sin(dip),u(1,:).*cos(strike)+u(2,:).*sin(strike));

    %Strikes are agreed to be positive
    %走向通常被认为是正值
    strikel0=find(strike<0);
    strike(strikel0)=strike(strikel0)+2*pi;
    fp=[strike' dip' rake'];

    %convert to Degree
    %将弧度转为角度
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