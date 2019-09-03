function newstruct = addfield(struct1,fields1,struct2,fields2)

struct1 = jaaba;
fields1 = {'roll_beg_short','bStart'};
len1 = length(fields1);
getfield(struct1,fields1{:})

temp = struct1;
for ii = 1:len1
    temp = temp.(fields1{ii});
end
