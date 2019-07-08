%demandtime = [0, 3600, 7200, 10800, 14400];
%demanddata = [1, .75, .5, 1, 1];

% test Visura
demandtime = [0, 3600, 7200, 10800 14400 18000];
demanddata = [1, 0.8, 0.5, 0.5, 0.5, 1];

demand = timeseries(demanddata,demandtime);

% area([0 3800], [L L])
% hold on; area((OTSG_mux.time)/1, (OTSG_mux.data(:,28) + OTSG_mux.data(:,30)))
% hold on; area((OTSG_mux.time)/1, OTSG_mux.data(:,30))
% plot(OTSG_mux.time - 1000, OTSG_mux.data(:,[25 28 30]))