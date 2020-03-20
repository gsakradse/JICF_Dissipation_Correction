function out = list_mats(pth)

% Function which takes a directory path and returns a cell array of the .mat
% files in that directory

out = dir(pth);
out = {out.name}';
out = out(contains(out, '.mat'));

end