% SPECSHIFT signal time shift by non-integer lag via spectral domain
%
% SYNTAX        datout = specshift(signal,shift)
%
% PURPOSE        Subsample time-shift, via phase shift in spectral domain 
%                        时移二次抽样，是在频率域做相移
%
% INPUT                signal        can be a matrix signal(time,chanel)
%                shift        in sample units, can be non-integer 
%                        can be a vector (one shift for each signal chanel)
%
% NOTE                In principle the input signal must be band-limited
%                to < Nyquist frequency (1/2 sampling frequency)
%                => preprocess by taper then low-pass filter
%                (this is usually done by the instrument anti-alias filter)

% Jean-Paul Ampuero        ampuero@erdw.ethz.ch

function datout = specshift(datin,shift)

if any(size(datin)==1), datin = datin(:); end
shift=shift(:);      % change shift to a column vector

nshift = length(shift);
[Nin,nsis] = size(datin);
if nshift>1 & nsis>1 & nshift~=nsis, error('arguments have wrong size'); end
% nextpow2(N) returns the first P such that 2.^P >= abs(N)
N = 2^(nextpow2( Nin +max(abs(shift)) )); % +? to avoid wraparound
fdatin = fft(datin,N);

ik = 2i*pi*(0:N/2).'/N;      % 0:N/2 is row vector, .' is common transposition, different from '
ik(N/2+2:N)=-ik(N/2:-1:2);

if nsis==1 & nshift==1   % only one data
  
  fshift = exp(ik*shift) ;
  fshift(N/2+1)=real(fshift(N/2+1));         % must be real at the Nyquist frequency (->output is real),
                                          % although irrelevant when the input has been low-pass filtered
  datout = real(ifft(fshift .*fdatin));
  datout = datout(1:Nin) ;
  
else
  
  fshift = exp(ik*shift') ;  % [N,ns]
  if nsis==1 & nshift>1, fdatin = repmat(fdatin,1,nshift); end
  if nshift==1 & nsis>1, fshift = repmat(fshift,1,nsis);   end
  fshift(N/2+1,:)=real(fshift(N/2+1,:));
  fshift = fshift.*fdatin;
  datout = real(ifft( fshift ));
  datout = datout(1:Nin,:);
  
end
