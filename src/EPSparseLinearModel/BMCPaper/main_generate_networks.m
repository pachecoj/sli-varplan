rng(0);

N = 50; dir = 'data/networks_SW'; model = 'smallworld';
if ~exist(dir,'dir'), mkdir(dir); end;

for r = 1:50
    fprintf('Network = %i\n',r);
    G = generate_network(N,model);
    save(sprintf('%s/N=%i_r=%i.mat',dir,N,r),'G')
end
