function s = rand_str(len)
s = [];
for i = 1:len
    s = [s char('a'+floor(rand*26))];
end