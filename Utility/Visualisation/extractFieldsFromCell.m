function data = extractFieldsFromCell(cell)
names = fields(cell{1});
mat = cell2mat(cell);
for k = 1:numel(names)
    data.(names{k}) = cell2mat({mat.(names{k})});
end
end