actions = {'prepare','run'};
settings = {'design','noise','inputscale','inputshape'};
for i = 2 % length(actions)
    for j = 1 %:length(settings)
        main(actions{i},settings{j});
    end
end
