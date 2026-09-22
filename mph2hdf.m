%MPH2HDF Extraction COMSOL -> HDF5 (LiveLink for MATLAB).
%
%   Pour chaque combinaison de paramètres : résolution du modèle COMSOL,
%   extraction des matrices réduites K, M et des positions des noeuds, puis
%   écriture dans le fichier HDF5 :
%
%   /<GROUP_NAME>/<id>      attributs : l, w, ..., dofs, nodes
%       K, M   (n_dofs, n_dofs)
%       XYZ    (n_nodes, 3)    positions en um, noeuds '1', '2', ...
%
%   Une configuration déjà présente (mêmes paramètres) est écrasée.
%   Formes indiquées telles que lues en Python (h5py).
%   À lancer depuis « COMSOL Multiphysics with MATLAB ».

H5_FILE     = 'aps.h5';                         % fichier de sortie
DATASET_ROM = 'rom1_n_rfc1_solid_1';            % solution du modèle réduit
SYS_MATRIX  = 'sys1';                           % noeud « Matrice système »

%% ──────────────────────────────── ZONE UTILISATEUR ──────────────────────
MPH_FILE   = 'C:\Users\GM287120\Desktop\Comsol\SuperElements\p110.mph';
GROUP_NAME = 'p110';

dofs = {'eta1', ...
        '1.ux', '1.uy', '1.uz', '1.rx', '1.ry', '1.rz', ...
        '2.ux', '2.uy', '2.uz', '2.rx', '2.ry', '2.rz'};

nodes = {'comp1.solid.att1', ...                % noeud 1
         'comp1.solid.att2'};                   % noeud 2

parameters = struct('l', [15 30 130], ...
                    'w', 3);

S = [1, repmat([1 1 1 0.5e6 0.5e6 0.5e6], 1, 2)];   % mise à l'échelle des DDL

%% ────────────────────────────────── EXTRACTION ──────────────────────────
assert(numel(S) == numel(dofs), 'S a %d termes pour %d DDL', numel(S), numel(dofs));

combos = balayage(parameters);
meta   = struct('dofs',  {dofs}, ...
                'nodes', {cellstr(string(1:numel(nodes)))});   % '1', '2', ...

disp('Chargement du modèle...')
model   = mphload(MPH_FILE);
studies = model.study.tags;
study   = model.study(studies(1));              % première étude du modèle

t0 = tic;
for i = 1:numel(combos)
    t = tic;
    p = combos(i);
    fprintf('[%d/%d] %s', i, numel(combos), jsonencode(p));

    % Résolution
    for f = fieldnames(p).'
        model.param.set(f{1}, num2str(p.(f{1}), 16));
    end
    study.run;

    % Extraction
    sys  = model.result.numerical(SYS_MATRIX);
    data = struct('K',   reducedMatrix(sys, DATASET_ROM, 'stiffness', S), ...
                  'M',   reducedMatrix(sys, DATASET_ROM, 'mass', S), ...
                  'XYZ', attCoords(model, nodes));

    % Écriture
    id = save2hdf(H5_FILE, GROUP_NAME, p, data, meta);
    fprintf(' -> %s/%s (%.1f s)\n', GROUP_NAME, id, toc(t));
end
fprintf('Terminé en %.1f s -> %s\n', toc(t0), H5_FILE);

%% ─────────────────────────────── FONCTIONS : CALCUL ─────────────────────
function combos = balayage(params)
%BALAYAGE Toutes les combinaisons des paramètres -> tableau de structs.
%   struct('l', [15 30], 'w', 3)  ->  (l=15, w=3), (l=30, w=3)
names = fieldnames(params);
vals  = zeros(1, 0);                            % une ligne par combinaison
for k = 1:numel(names)
    v    = params.(names{k})(:);
    % chaque combinaison existante est déclinée pour chaque valeur de v
    vals = [repelem(vals, numel(v), 1), repmat(v, size(vals, 1), 1)];
end
combos = cell2struct(num2cell(vals), names, 2);
end

function A = reducedMatrix(sys, dset, kind, S)
%REDUCEDMATRIX Matrice réduite du ROM ('stiffness' ou 'mass'), mise à l'échelle par S.
sys.set('solution', dset);
sys.set('reducedmodelmatrix', kind);
sys.set('format', 'filled');
A = sys.getReal();
n = numel(S);
assert(isequal(size(A), [n n]), 'COMSOL renvoie une matrice %dx%d pour %d DDL', size(A), n);
A = cleanMatrix(S(:) .* A .* S(:).');           % A(i,j) * S(i) * S(j)
end

function A = cleanMatrix(A)
%CLEANMATRIX Symétrise A et annule les termes hors diagonale négligeables,
%   c.-à-d. |A(i,j)| < 1e-8 * sqrt(|A(i,i) * A(j,j)|).
A     = (A + A.') / 2;
d     = abs(diag(A));
ref   = sqrt(d * d.');
noise = abs(A) < 1e-8 * ref | ref == 0;
noise(logical(eye(size(A)))) = false;           % la diagonale est conservée
A(noise) = 0;
end

function xyz = attCoords(model, tags)
%ATTCOORDS Positions (um) des attachements : une ligne [x y z] par noeud.
xyz = zeros(numel(tags), 3);
for i = 1:numel(tags)
    expr = strcat(tags{i}, {'.xcx', '.xcy', '.xcz'});
    [x, y, z] = mphglobal(model, expr, 'unit', {'um', 'um', 'um'});
    xyz(i, :) = [x(1), y(1), z(1)];
end
xyz(abs(xyz) < 1e-15) = 0;                      % zéros numériques
end

%% ──────────────────────────────── FONCTIONS : HDF5 ──────────────────────
function id = save2hdf(h5File, group, p, data, meta)
%SAVE2HDF Écrit une configuration dans /<group>/<id> et renvoie son id.
grp = ['/' group];
id  = findConfig(h5File, grp, p);
cfg = [grp '/' id];
for f = fieldnames(data).'                      % datasets K, M, XYZ
    writeDataset(h5File, [cfg '/' f{1}], data.(f{1}));
end
for f = fieldnames(p).'                         % attributs l, w, ...
    h5writeatt(h5File, cfg, f{1}, p.(f{1}));
end
for f = fieldnames(meta).'                      % attributs dofs, nodes
    writeStrAttr(h5File, cfg, f{1}, meta.(f{1}));
end
end

function id = findConfig(h5File, grp, p)
%FINDCONFIG Id de la configuration de mêmes paramètres, sinon premier id libre.
try
    info = h5info(h5File, grp);
    cfgs = info.Groups;
catch
    cfgs = [];                                  % fichier ou groupe pas encore créé
end
ids = 0;
for k = 1:numel(cfgs)
    id = cfgs(k).Name(numel(grp) + 2:end);      % '/p110/12' -> '12'
    if sameParams(cfgs(k).Attributes, p), return, end
    ids(end + 1) = str2double(id); %#ok<AGROW>
end
id = num2str(max(ids) + 1);
end

function tf = sameParams(attrs, p)
%SAMEPARAMS Vrai si les attributs du groupe portent les valeurs de p (comme numpy.isclose).
isclose = @(a, b) abs(a - b) <= 1e-8 + 1e-5 * abs(b);
tf = false;
if isempty(attrs), return, end                  % groupe incomplet
for f = fieldnames(p).'
    k = strcmp({attrs.Name}, f{1});
    if ~any(k) || ~isclose(attrs(k).Value, p.(f{1})), return, end
end
tf = true;
end

function writeDataset(h5File, ds, x)
%WRITEDATASET Écrit la matrice x (double, gzip), en écrasant le dataset existant.
%   x est transposé : MATLAB range par colonnes, h5py lit par lignes.
x = x.';
try
    h5create(h5File, ds, size(x), 'ChunkSize', size(x), 'Deflate', 4);
catch err                                       % déjà créé : on réécrit dedans
    if ~strcmp(err.identifier, 'MATLAB:imagesci:h5create:datasetAlreadyExists'), rethrow(err), end
end
h5write(h5File, ds, x);
end

function writeStrAttr(h5File, loc, name, str)
%WRITESTRATTR Attribut « liste de textes » au format h5py (UTF-8, longueur variable).
%   h5writeatt ne sait pas écrire ce type sur toutes les versions de MATLAB.
fid   = H5F.open(h5File, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid   = H5G.open(fid, loc);
type  = H5T.copy('H5T_C_S1');                   % chaîne de caractères...
H5T.set_size(type, 'H5T_VARIABLE');             % ...de longueur variable
H5T.set_cset(type, H5ML.get_constant_value('H5T_CSET_UTF8'));
space = H5S.create_simple(1, numel(str), []);   % tableau 1D
try, H5A.delete(gid, name); catch, end          % supprime l'ancien s'il existe
attr  = H5A.create(gid, name, type, space, 'H5P_DEFAULT');
H5A.write(attr, type, str);
H5A.close(attr); H5S.close(space); H5T.close(type); H5G.close(gid); H5F.close(fid);
end
