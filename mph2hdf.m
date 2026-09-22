%MPH2HDF Extraction COMSOL -> HDF5 (LiveLink for MATLAB).
%
%   /<GROUP_NAME>        attributs : dofs (n_dofs,), nodes (n_nodes,)
%       /<config_id>     attributs : l, w, ...
%           K, M   (n_dofs, n_dofs)
%           XYZ    (n_nodes, 3)  positions en um, ordre de `nodes`
%
%   Les positions sont lues sur les attachements (att.xcx/.xcy/.xcz) après
%   résolution : elles suivent le balayage paramétrique.
%
%   Les formes ci-dessus sont celles lues en Python (h5py, ordre C) : les
%   tableaux sont transposés à l'écriture pour produire le même fichier que
%   mph2hdf.py.
%
%   À lancer depuis « COMSOL Multiphysics with MATLAB » (ou après mphstart).

H5_FILE     = 'aps.h5';
DATASET_ROM = 'rom1_n_rfc1_solid_1';
SYS_MATRIX  = 'sys1';

%% ──────────────────────────────── ZONE UTILISATEUR ──────────────────────
MPH_FILE   = 'C:\Users\GM287120\Desktop\Comsol\SuperElements\p110.mph';
GROUP_NAME = 'p110';

dofs = {'eta1', ...
        'n1.ux', 'n1.uy', 'n1.uz', 'n1.rx', 'n1.ry', 'n1.rz', ...
        'n2.ux', 'n2.uy', 'n2.uz', 'n2.rx', 'n2.ry', 'n2.rz'};

nodes = struct('n1', 'comp1.solid.att1', ...
               'n2', 'comp1.solid.att2');

parameters = struct('l', [15 30 130], ...
                    'w', 3);

S = [1, repmat([1 1 1 0.5e6 0.5e6 0.5e6], 1, 2)];
%% ──────────────────────────────────── BACKEND ───────────────────────────
assert(numel(S) == numel(dofs), 'S a %d termes pour %d DDL', numel(S), numel(dofs));

combos   = balayage(parameters);
grp      = ['/' GROUP_NAME];
grpAttrs = struct('dofs', {dofs}, 'nodes', {fieldnames(nodes)});

disp('Chargement du modèle COMSOL...')
model   = mphload(MPH_FILE);
studies = model.study.tags;
study   = model.study(studies(1));

fprintf('\nBalayage paramétrique (%d configurations) -> %s\n', numel(combos), H5_FILE);

for i = 1:numel(combos)
    p = combos(i);
    fprintf('\n[%d/%d] %s\n', i, numel(combos), jsonencode(p));

    for name = fieldnames(p).'
        model.param.set(name{1}, num2str(p.(name{1}), 16));
    end
    study.run;

    sys  = model.result.numerical(SYS_MATRIX);
    XYZ  = attCoords(model, struct2cell(nodes));
    data = struct('K',   reducedMatrix(sys, DATASET_ROM, 'stiffness', S), ...
                  'M',   reducedMatrix(sys, DATASET_ROM, 'mass', S), ...
                  'XYZ', XYZ);

    fprintf('    XYZ (um) : %s\n', mat2str(round(XYZ, 3)));
    fprintf('    -> %s/%s\n', GROUP_NAME, save2hdf(H5_FILE, grp, grpAttrs, p, data));
end

com.comsol.model.util.ModelUtil.remove(model.tag);
disp('Extraction terminée !')

%% ──────────────────────────────── FONCTIONS ─────────────────────────────
function combos = balayage(params)
%BALAYAGE Produit cartésien des paramètres -> tableau de structs {nom: valeur}.
names = fieldnames(params);
vals  = zeros(1, 0);                            % une combinaison vide = modèle nominal
for k = 1:numel(names)
    v    = params.(names{k})(:);
    vals = [repelem(vals, numel(v), 1), repmat(v, size(vals, 1), 1)];
end
combos = cell2struct(num2cell(vals), names, 2);
end

function A = cleanMatrix(A, tol)
%CLEANMATRIX Symétrise et supprime le bruit numérique (diagonale conservée).
if nargin < 2, tol = 1e-8; end
A     = (A + A.') / 2;
d     = abs(diag(A));
denom = sqrt(d * d.');
A((abs(A) < tol * denom | denom == 0) & ~eye(size(A))) = 0;
end

function xyz = attCoords(model, tags)
%ATTCOORDS Positions (um) des attachements -> matrice (n_nodes, 3).
expr = cellstr(string(tags(:)).' + [".xcx"; ".xcy"; ".xcz"]);   % 3 x n_nodes
vals = cell(size(expr));
[vals{:}] = mphglobal(model, expr(:).', 'unit', repmat({'um'}, 1, numel(expr)));
xyz = cellfun(@(v) v(1), vals).';
xyz(abs(xyz) < 1e-15) = 0;
end

function A = reducedMatrix(sys, dset, kind, S)
%REDUCEDMATRIX Matrice réduite du ROM ('stiffness' ou 'mass'), mise à l'échelle par S.
sys.set('solution', dset);
sys.set('reducedmodelmatrix', kind);
sys.set('format', 'filled');
A = sys.getReal();
n = numel(S);
assert(isequal(size(A), [n n]), 'COMSOL renvoie une matrice %dx%d pour %d DDL', size(A), n);
A = cleanMatrix(S(:) .* A .* S(:).');
end

function id = save2hdf(h5File, grp, grpAttrs, p, data)
%SAVE2HDF Écrit une configuration dans le HDF5, renvoie son identifiant.
[id, isNew] = configId(h5File, grp, p);
cfg = [grp '/' id];
for name = fieldnames(data).'
    writeDataset(h5File, [cfg '/' name{1}], data.(name{1}));
end
if isNew
    for name = fieldnames(p).'
        h5writeatt(h5File, cfg, name{1}, p.(name{1}));
    end
end
for name = fieldnames(grpAttrs).'
    writeStrAttr(h5File, grp, name{1}, grpAttrs.(name{1}));
end
end

function [id, isNew] = configId(h5File, grp, p)
%CONFIGID Configuration existante de mêmes paramètres, sinon nouvel identifiant.
try
    info = h5info(h5File, grp);
    cfgs = info.Groups;
catch
    cfgs = [];                                  % fichier ou groupe inexistant
end
isNew = false;
ids   = 0;
for c = reshape(cfgs, 1, [])
    id = c.Name(numel(grp) + 2:end);            % '/p110/12' -> '12'
    if ~isempty(regexp(id, '^\d+$', 'once'))
        if sameParams(c, p), return, end
        ids(end + 1) = str2double(id); %#ok<AGROW>
    end
end
id    = num2str(max(ids) + 1);
isNew = true;
end

function tf = sameParams(info, p)
%SAMEPARAMS Vrai si le groupe (h5info) porte les valeurs de p, à np.isclose près.
isclose = @(a, b) abs(a - b) <= 1e-8 + 1e-5 * abs(b);
a = info.Attributes;
if isempty(a), a = struct('Name', {}, 'Value', {}); end
tf = true;
for name = fieldnames(p).'
    k  = strcmp({a.Name}, name{1});
    tf = tf && any(k) && isclose(a(k).Value, p.(name{1}));
end
end

function writeDataset(h5File, ds, x)
%WRITEDATASET Écrit x en double compressé (gzip), en écrasant le dataset existant.
x = x.';                                        % ordre colonnes -> ordre C (h5py)
try
    h5create(h5File, ds, size(x), 'ChunkSize', size(x), 'Deflate', 4);
catch err
    if ~strcmp(err.identifier, 'MATLAB:imagesci:h5create:datasetAlreadyExists'), rethrow(err), end
end
h5write(h5File, ds, x);
end

function writeStrAttr(h5File, loc, name, str)
%WRITESTRATTR Attribut tableau de chaînes UTF-8 de longueur variable (comme h5py).
fid   = H5F.open(h5File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid   = H5G.open(fid, loc);
type  = H5T.copy('H5T_C_S1');
H5T.set_size(type, 'H5T_VARIABLE');
H5T.set_cset(type, H5ML.get_constant_value('H5T_CSET_UTF8'));
space = H5S.create_simple(1, numel(str), []);
try, H5A.delete(gid, name); catch, end          % remplace l'attribut s'il existe
attr  = H5A.create(gid, name, type, space, 'H5P_DEFAULT');
H5A.write(attr, type, str);
H5A.close(attr); H5S.close(space); H5T.close(type); H5G.close(gid); H5F.close(fid);
end
