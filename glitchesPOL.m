function nzeros=glitchesPOL(trace,nwin,winlen,winoff,igstart,sps,STAsoff)
% In function glitchesPOL, we're counting CONSECUTIVE zeros.  @100sps.

nzeros=zeros(nwin,1);
spsfactor=100/sps; %assumes that POLARIS stations record at 100 sps, and that igstart is given in sps samples.
igstart=round((igstart-1)*spsfactor)+1+STAsoff;
a=STAsoff
winoff=round(winoff*spsfactor);
winlen=round(winlen*spsfactor);
% Could expand window, in case filtered glitch extends into this window from adjacent zeros.
for n=1:nwin
    istart=igstart+(n-1)*winoff;
    iend=istart+winlen-1;
    nzeros(n)=0;
    ii=istart-1;
    while ii<iend-1
        ii=ii+1;
        nz=0;
        while (abs(trace(ii)) <= 1.e-07 && ii<iend-1)
            nz=nz+1;
            ii=ii+1;
        end
        nzeros(n)=max(nzeros(n),nz); 
%         if nzeros(n) >= 5
%             a=[n nzeros(n) ii]
%             b=trace(istart:iend)'
%         end
    end
end

%     if STAsoff > -1
%         STAscrot(1:tracelen-STAsoff)=STAscrot(STAsoff+1:tracelen);
%         STAscrot(tracelen-STAsoff+1:tracelen)=0;
%     else
%         STAscrot(-STAsoff+1:tracelen)=STAscrot(1:tracelen+STAsoff);
%         STAscrot(1:-STAsoff)=0;
%     end
