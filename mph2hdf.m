%MPH2HDF Extraction COMSOL -> HDF5 (à lancer depuis « COMSOL with MATLAB »).
%   Pour chaque combinaison de paramètres, écrit dans H5_FILE le groupe
%   /GROUP_NAME/<id> : attributs l, w, ..., dofs, datasets K, M (n_dofs, n_dofs)
%   et XYZ (n_nodes, 3) en um, noeuds dans l'ordre de `nodes`.
%   Une configuration déjà présente (mêmes paramètres) est remplacée.

H5_FILE     = 'aps.h5';
DATASET_ROM = 'rom1_n_rfc1_solid_1';
SYS_MATRIX  = 'sys1';

%% Zone utilisateur
MPH_FILE   = 'C:\Users\GM287120\Desktop\Comsol\SuperElements\p110.mph';
GROUP_NAME = 'p110';

dofs = {'eta1', ...
        '1.ux', '1.uy', '1.uz', '1.rx', '1.ry', '1.rz', ...
        '2.ux', '2.uy', '2.uz', '2.rx', '2.ry', '2.rz'};

nodes = {'comp1.solid.att1', ...                % noeud 1
         'comp1.solid.att2'};                   % noeud 2

parameters = struct('l', [15 30 130], ...
                    'w', 3);

S = [1, repmat([1 1 1 0.5e6 0.5e6 0.5e6], 1, 2)];

%% Extraction
assert(numel(S) == numel(dofs), 'S a %d termes pour %d DDL', numel(S), numel(dofs));

names  = fieldnames(parameters);
grids  = struct2cell(parameters);
[grids{:}] = ndgrid(grids{:});                  % toutes les combinaisons
combos = cell2mat(cellfun(@(g) g(:), grids', 'UniformOutput', false));   % une ligne par combinaison

model   = mphload(MPH_FILE);
studies = model.study.tags;

t0 = tic;
for i = 1:size(combos, 1)
    t = tic;
    fprintf('[%d/%d] %s', i, size(combos, 1), strjoin(string(names') + "=" + combos(i, :), ', '));

    for k = 1:numel(names)
        model.param.set(names{k}, num2str(combos(i, k), 16));
    end
    model.study(studies(1)).run;

    sys  = model.result.numerical(SYS_MATRIX);
    data = struct('K',   reducedMatrix(sys, DATASET_ROM, 'stiffness', S), ...
                  'M',   reducedMatrix(sys, DATASET_ROM, 'mass', S), ...
                  'XYZ', attCoords(model, nodes));

    id = save2hdf(H5_FILE, GROUP_NAME, names, combos(i, :), data, dofs);
    fprintf(' -> %s/%d (%.1f s)\n', GROUP_NAME, id, toc(t));
end
fprintf('Terminé en %.1f s\n', toc(t0));

%% Fonctions
function A = reducedMatrix(sys, dset, kind, S)
%REDUCEDMATRIX Matrice réduite du ROM ('stiffness' ou 'mass') mise à l'échelle par S,
%   symétrisée, termes hors diagonale < 1e-8 * sqrt(|A(i,i) * A(j,j)|) mis à zéro.
sys.set('solution', dset);
sys.set('reducedmodelmatrix', kind);
sys.set('format', 'filled');
A   = S(:) .* sys.getReal() .* S(:).';
A   = (A + A.') / 2;
ref = sqrt(abs(diag(A)) * abs(diag(A)).');
A((abs(A) < 1e-8 * ref | ref == 0) & ~eye(size(A))) = 0;
end

function xyz = attCoords(model, tags)
%ATTCOORDS Positions (um) des attachements : une ligne [x y z] par noeud.
xyz = zeros(numel(tags), 3);
for i = 1:numel(tags)
    [x, y, z] = mphglobal(model, strcat(tags{i}, {'.xcx', '.xcy', '.xcz'}), 'unit', {'um', 'um', 'um'});
    xyz(i, :) = [x(1), y(1), z(1)];
end
xyz(abs(xyz) < 1e-15) = 0;
end

function id = save2hdf(h5File, group, names, vals, data, dofs)
%SAVE2HDF Écrit une configuration dans /group/<id> et renvoie id.
[id, found] = findConfig(h5File, group, names, vals);
cfg = sprintf('/%s/%d', group, id);
for f = fieldnames(data)'                       % transposé : h5py lit par lignes
    x = data.(f{1}).';
    if ~found, h5create(h5File, [cfg '/' f{1}], size(x)); end
    h5write(h5File, [cfg '/' f{1}], x);
end
for k = 1:numel(names)
    h5writeatt(h5File, cfg, names{k}, vals(k));
end
writeStrAttr(h5File, cfg, 'dofs', dofs);
end

function [id, found] = findConfig(h5File, group, names, vals)
%FINDCONFIG Id de la configuration de mêmes paramètres, sinon premier id libre.
try, info = h5info(h5File, ['/' group]); catch, info.Groups = []; end
id = 0;
for g = info.Groups'
    n = str2double(g.Name(numel(group) + 3:end));   % '/p110/12' -> 12
    [ok, k] = ismember(names, {g.Attributes.Name});
    found = all(ok) && all(abs([g.Attributes(k).Value] - vals) <= 1e-8 + 1e-5 * abs(vals));
    if found, id = n; return, end
    id = max(id, n);
end
[id, found] = deal(id + 1, false);
end

function writeStrAttr(h5File, loc, name, str)
%WRITESTRATTR Écrit l'attribut « liste de textes » str (format h5py).
fid = H5F.open(h5File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid = H5G.open(fid, loc);
typ = H5T.copy('H5T_C_S1');
H5T.set_size(typ, 'H5T_VARIABLE');
spc = H5S.create_simple(1, numel(str), []);
try, H5A.delete(gid, name); catch, end          % remplace l'attribut existant
att = H5A.create(gid, name, typ, spc, 'H5P_DEFAULT');
H5A.write(att, typ, str);
H5A.close(att); H5S.close(spc); H5T.close(typ); H5G.close(gid); H5F.close(fid);
end
