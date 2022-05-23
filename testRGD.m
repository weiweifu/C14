fname = '/global/cfs/cdirs/e3sm/weiweif/MPAS_IRF/TestCases_IFVH/EC30to60E2r2/timeAverage.mon.00.nc';

vn = 'timeMonthly_avg_debugTracers_tracer1';
vn = 'timeMonthly_avg_activeTracers_temperature';

tracer1 = ncread(fname,vn,[1 1 1],[inf inf inf]);
sst = tracer1(1,:);

dst_field = repmat(0.0,rgd.n_b,1);
% Apply weights
for i=1:rgd.n_s
  % dst_field(:,rgd.row(i)) = rgd.S(i)*tracer1(:,rgd.col(i));
  dst_field(rgd.row(i)) = dst_field(rgd.row(i)) + rgd.S(i).*sst(rgd.col(i));
end

% Adjust destination field by fraction
for i=1:rgd.n_b
  if rgd.frac_b(i) > 0.0 
     dst_field(i) = dst_field(i)./rgd.frac_b(i);
  end
end

S = sparse(double(rgd.row(:)),double(rgd.col(:)),rgd.S(:),double(rgd.n_b),double(rgd.n_a));

frac_b = rgd.frac_b;
frac_b(frac_b <= 0) = nan;

sst2 = S*double(sst')./frac_b;
