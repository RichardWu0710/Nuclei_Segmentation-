function Combo_Output = Combo_Cell(old,new)

old = old(1:end-1,1:end-1);
new = new(1:end-1,1:end-1);

Genestr = [old(2:end,1); new(2:end,1)];
Gene = unique(Genestr);

oldrow = old(1,1:end);
oldcol = old(1:end,1);
newrow = new(1,1:end);
newcol = new(1:end,1);

Cell = [old(1,2:end) new(1,2:end)];

ComboColLen = length(Gene)+2;
ComboRowLen = length(Cell)+2;

Combo = cell(ComboColLen,ComboRowLen);
Combo{1,1} = 'Gene';
Combo{end,1} = 'Total';
Combo{1,end} = 'Total';

for i = 2:ComboColLen-1
    Combo{i,1} = Gene{i-1};
end

for i = 2:ComboRowLen-1
    Combo{1,i} = Cell{i-1};
end

Comborow = Combo(1,1:end);
Combocol = Combo(1:end,1);
for i = 1:ComboRowLen
    for j = 1:ComboColLen
        if isempty(Combo{j,i})
            Combo{j,i} = 0;
        end
    end
end

for i = 2:length(oldrow)
   for j = 2:length(oldcol)
       singene = old{j,1};
       yind = strcmp(Combocol,singene);
       yind = find(yind);
       Combo{yind,i} = Combo{yind,i}+old{j,i};
   end
end

for i = 2:length(newrow)
   for j = 2:length(newcol)
       singene = new(j,1);
       yind = strcmp(Combocol,singene);
       yind = find(yind);
       Combo{yind,i-1+length(oldrow)} = Combo{yind,i-1+length(oldrow)}+new{j,i};
   end
end

% Sums of genes for each cell
for i = 2:ComboRowLen
    sum = 0;
    for j = 2:ComboColLen
       sum = sum + Combo{j,i};
    end
    Combo{end,i} = sum;
end

% Total number of each type of genes
for i = 2:ComboColLen
    sum = 0;
    for j = 2:ComboRowLen
       sum = sum + Combo{i,j};
    end
    Combo{i,end} = sum;
end


Combo_Output = Combo;
end