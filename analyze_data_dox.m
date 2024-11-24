
% typical script to analyze a set of data

clc;
%% image paths

p = '/Users/aliarra/Documents/WDoxExp/';

% For DOX
name = {'WT_glyTP2_LD_01_ISeg.tif',...
    'WT_glyTP2_LD_02_ISeg.tif','WT_glyTP2_LD_03_ISeg.tif'};

namer = {'WT_glyTP2_LD_01_ISeg.tif',...
    'WT_glyTP2_LD_02_ISeg.tif','WT_glyTP2_LD_03_ISeg.tif'};

cellsname = {'WT_glyTP2_cells_01_er.tif',...
    'WT_glyTP2_cells_02_er.tif','WT_glyTP2_cells_03_er.tif'};


%%
tic;
for i=1:3
    data(i) = dataclass([p,name{i}]);
    data(i) = data(i).loadOrgImage();
    data(i).CellPath = [p,cellsname{i}];
    data(i) = data(i).loadCellImage();
    data(i).AuxPath = [p,namer{i}];
    data(i) = data(i).loadAuxImage();
end
toc

%% params 

SIGMA = 0.5;
THRESHOLD = 0.1;

MINSIZE = 3; MAXSIZE = 10000; %3 for globular & 6 for non
MININT = 0;
MINCELLSIZE = 100; 

%% analyze
tic;
for i=1:3
    statistics(i) = analyze(data(i).OrgImage,data(i).AuxImage,data(i).CellImage,...
        SIGMA,THRESHOLD,MINSIZE,MAXSIZE,MININT,MINCELLSIZE);
end
toc

statistics_merged = statistics(1);
% concatenate into data
statistics_merged = CatStructFields(statistics(1),statistics(2),1);
for i=1:2
    statistics_merged = CatStructFields(statistics_merged,statistics(i),1);
end
% 
statistics_merged = CatStructFields(statistics_merged,statistics(3),1);
for i=2:2
    statistics_merged = CatStructFields(statistics_merged,statistics(i),1);
end

%% function
function data = analyze(orgs,orgs_raw,cells,SIGMA,THRESHOLD,MINSIZE,MAXSIZE,MININT,MINCELLSIZE)

    numz = size(orgs,3);
    px_c = .433332;

    %% binarize cells
    bw = cells > 0;

cells_bwinit = imbinarize(cells);

D = -bwdist(~cells_bwinit);
Ld = watershed(D);
bw2 = cells_bwinit;
bw2(Ld == 0) = 0;
mask = imextendedmin(D,2);
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = cells_bwinit;
bw3(Ld2 == 0) = 0;
bw=bw3;

fileID=fopen('transform_inverse.txt','r');  %read in transform matrix
formatSpec = '%f';
sizeA = [3 3];
A = fscanf(fileID,formatSpec,sizeA);
A = inv(A);
A(7)=0;
A(8)=0;
A(9)=1;
tform = affine2d(A);                        %define affine matrix

outView = imref2d([512 512]);       %set output ref for affine transform
mapoforgsandcells = labeloverlay(orgs(:,:,18),bw);


    %% segment cells
    [cell_bounds,L] = bwboundaries(bw,4,'noholes'); % why 4? chosen from experience (ie staring at images)
    NUM_CELLS = length(cell_bounds);   % number of cells in population; includes cells that may later be filtered out

    cell_sizes = zeros([length(cell_bounds),1]);  % we hold our cell sizes here
    cell_coords = cell([length(cell_bounds),1]);   % we hold our cell boundary coords here (rectangle)
    cell_volumes = zeros([length(cell_bounds),1]);

    % iter thru cells 
    for c=1:NUM_CELLS
        % grab cell boundary coordinates:
        these_bounds = cell_bounds{c};
        % draw rectangle around cell:
        xmin = min(these_bounds(:,2))-1; xmax = max(these_bounds(:,2))+1;
        ymin = min(these_bounds(:,1))-1; ymax = max(these_bounds(:,1))+1;
        % store rectangle:
        cell_coords{c} = [xmin xmax; ymin ymax];

        if min(cell_coords{c}(:)) < 1 || max(cell_coords{c}(:)) > size(bw,1)   % image borders are unreliable
            cell_sizes(c) = NaN;
            continue
        else
            this_c = bw(ymin:ymax,xmin:xmax,:);
            if sum(this_c(:)) < MINCELLSIZE
                cell_sizes(c) = NaN;
                continue
            else
                cell_sizes(c) = sum(this_c(:));
            end
        end
cell_label = bwlabel(L,4);
stat=regionprops(L,'centroid','Orientation','MajorAxisLength','MinorAxisLength','area','boundingbox','eccentricity', 'pixellist');
stats = regionprops('table',L,'centroid','MajorAxisLength','MinorAxisLength','area','boundingbox','eccentricity','pixellist');
cell_area = cat(1,stat.Area);
centroids = cat(1,stat.Centroid);
major = cat(1,stat.MajorAxisLength);
minor = cat(1,stat.MinorAxisLength);

cell_sizeum = cell_sizes.*(.414320)^2;
height = major.*minor;
height = sqrt(height);
cell_volum = cell_sizeum.*height*.2;

cell_volumes = (4*pi/3).*(stats.MajorAxisLength./2).*(stats.MinorAxisLength./2).^2;
vol = cell_volumes.*cell_sizes;
cell_volumes= vol./cell_sizes;


    end
%  
% figure
% hold on
% imshow(L);
% plot(centroids(:,1),centroids(:,2),'b*')
% hold off

% figure
% hold on
% imshow(mapoforgsandcells);
% for kstat = 1:numel(stat);
%     cstat = stat(kstat).Centroid;
%     text(cstat(1), cstat(2), sprintf('%d', kstat), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'middle', 'Color', 'black', 'fontsize',12);
% end
% hold off
    %% partition cells into blocks
    seg_orgs = cell([NUM_CELLS,1]);
    seg_raw_orgs = cell([NUM_CELLS,1]);
    for c=1:NUM_CELLS
        if isnan(cell_sizes(c))
            continue
        else
            xmin = cell_coords{c}(1,1); xmax = cell_coords{c}(1,2);
            ymin = cell_coords{c}(2,1); ymax = cell_coords{c}(2,2);

            if xmax - xmin < 2 || ymax - ymin < 2
                continue
            else
                seg_orgs{c}(:,:,:) = orgs(ymin:ymax,xmin:xmax,:);
                seg_raw_orgs{c}(:,:,:) = orgs_raw(ymin:ymax,xmin:xmax,:);
            end
        end
    end

    %% now go back and measure their sizes
    seg_binary_orgs = cell([NUM_CELLS,1]);
    org_volumes = {};
    org_volumesum = {};
    org_ints = {};
    for c=1:NUM_CELLS
        if isempty(seg_orgs{c})
            org_volumes{c,1} = NaN;
            org_volumesum{c,1} = NaN;
            org_ints{c,1} = NaN;
            continue
        else
            this_cell = seg_orgs{c};
            int_mask = this_cell > MININT;
            this_raw_cell = seg_raw_orgs{c};
            this_bw = zeros(size(this_cell));
            pixels = ones([numz, 1]);
            pixelstds = ones([numz,1]);
            for z=1:numz
                x = this_cell(:,:,z);
                this_max = max(x(:));
                pixelstds(z) = std(x(:));
                if this_max < MININT
                    pixels(z) = 0;
                end
            end

            this_bw = (mylaplace(this_cell,SIGMA) > THRESHOLD) .* int_mask;
            for z=1:size(this_bw,3)
                if pixels(z) == 0
                    this_bw(:,:,z) = 0;
                end
            end
         
           seg_binary_orgs{c} = bwareaopen(this_bw,MINSIZE);
       

            this_cc = bwconncomp(this_bw,6); %6 for globular; 18 for large organelles
            for n=1:this_cc.NumObjects
                % grab an organelle from this cell:
                this_org = length(this_cc.PixelIdxList{n});
                % enforce size range on organelle 
                if this_org > MINSIZE && this_org < MAXSIZE
                    org_volumes{c,n} = this_org;
                    org_volumesum{c,n} = this_org.*0.0343;
                    org_ints{c,n} = sum(this_raw_cell(this_cc.PixelIdxList{n}));
                else
                    org_volumes{c,n} = NaN;
                    org_volumesum{c,n} = NaN;
                    org_ints{c,n} = NaN;
                end
            end
        end
    end

    for i=1:length(seg_binary_orgs);    
        bwconnstat(i) = bwconncomp(seg_binary_orgs{i},6);
            stats3D = regionprops3(bwconnstat(i),'Volume','SurfaceArea','Centroid');
            org_volcheck{i} = stats3D.Volume;
            org_SA{i} = stats3D.SurfaceArea;
            org_centroid{i} = stats3D.Centroid;
     end

    %% clean & save
    [rows,cols] = size(org_volumes);
    for i=1:rows
        for j=1:cols
            if isempty(org_volumes{i,j})
                org_volumes{i,j} = NaN;
                org_volumesum{i,j} = NaN;
                org_ints{i,j}=NaN;
            elseif org_volumes{i,j} < MINSIZE
                org_volumes{i,j} = NaN;
                org_volumesum{i,j} = NaN;
                org_ints{i,j} = NaN;
            end
        end
    end


    volumes = cell2mat(org_volumes); [M,N] = size(volumes);
    [~,col1] = sort(~isnan(volumes),2,'descend'); row1 = repmat(1:M,N,1)';
    restructured_indices = sub2ind(size(volumes),row1(:),col1(:));
    volumes = reshape(volumes(restructured_indices),M,N);
    volumes(:,all(isnan(volumes),1)) = [];
    
    ints = cell2mat(org_ints); [M,N] = size(ints);
    [~,col1] = sort(~isnan(ints),2,'descend'); row1 = repmat(1:M,N,1)';
    restructured_indices = sub2ind(size(ints),row1(:),col1(:));
    ints = reshape(ints(restructured_indices),M,N);
    ints(:,all(isnan(ints),1)) = [];

    volumesum = cell2mat(org_volumesum); [M,N] = size(volumesum);
    [~,col1] = sort(~isnan(volumesum),2,'descend'); row1 = repmat(1:M,N,1)';
    restructured_indices = sub2ind(size(volumesum),row1(:),col1(:));
    volumesum = reshape(volumesum(restructured_indices),M,N);
    volumesum(:,all(isnan(volumesum),1)) = [];

    
brows = cellfun(@numel,org_volcheck);
bcols = size(org_volcheck,2);
volcheck = zeros(max(brows),bcols);
for k =1:bcols
    volcheck(1:brows(k),k) = org_volcheck{k};
     surfacearea(1:brows(k),k) = org_SA{k};
end

volcheck = volcheck';
surfacearea = surfacearea';
surfacearea(volcheck < MINSIZE+1) = NaN;
volcheck(volcheck < MINSIZE+1) = NaN;

[M,N] = size(volcheck);
    [~,col1] = sort(~isnan(volcheck),2,'descend'); row1 = repmat(1:M,N,1)';
    restructured_indices = sub2ind(size(volcheck),row1(:),col1(:));
    volcheck = reshape(volcheck(restructured_indices),M,N);
    volcheck(:,all(isnan(volcheck),1)) = [];

    surfacearea = reshape(surfacearea(restructured_indices),M,N);
    surfacearea(:,all(isnan(surfacearea),1)) = [];
 
sphericity = ((pi)^(1/3)*(6.*volcheck).^(2/3))./surfacearea;

volfrac = sum(volumesum,2,'omitnan');
volfrac = volfrac./cell_volum;
    

    data.volumes = volumes;
    data.volumesum = volumesum; %org volumes in microns
    data.ints = ints;
    data.SA = surfacearea;
    data.volcheck = volcheck;
    data.cell_sizes = cell_sizes;
    data.cell_area = cell_area;
    PARAMS.THRESHOLD = THRESHOLD;
    PARAMS.SIGMA = SIGMA;
    PARAMS.MINSIZE = MINSIZE;
    PARAMS.MININT = MININT;
    data.seg_orgs = seg_orgs;
    data.seg_binary_orgs = seg_binary_orgs;
    data.THRESHOLD = THRESHOLD;
    data.SIGMA = SIGMA;
    data.PARAMS = PARAMS;
    data.bw = bw;
    out = volumes(any(~isnan(volumes),2),:); 
    PARAMS.N = size(out,1);
    data.cell_volumes = cell_volumes; %ellipsoid est
    data.cell_sizeum = cell_sizeum; %2D area of cells in microns
    data.cell_volum = cell_volum; %3D volume of cells in microns with cylinder est
    data.volfrac = volfrac; %volume fraction with micron data
    data.sphericity = sphericity;
    data.centroids = vertcat(centroids);  

end
