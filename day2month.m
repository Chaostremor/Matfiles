%nzeros=glitches2(PGCN,nwin,winlen,winoff,igstart,spsfactor);
function MO=day2month(jday,year)

daysinmonth=[31 28 31 30 31 30 31 31 30 31 30 31];
if abs(year-2004)<=1.e-6
    daysinmonth=[31 29 31 30 31 30 31 31 30 31 30 31];
end
MONTHS=char('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC');
MO=MONTHS(find(jday<=cumsum(daysinmonth),1),:);
