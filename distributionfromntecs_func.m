function [even_along_vess, random_ntecs, actual_ntecs, results_table] = distributionfromntecs_func(vessels_labeled_segs, binary_hs_image, skel_img, dapi_img, px_per_um, save_dir,sample_name)
tic
shortfile = sample_name;
display(['Analyzing nanoparticle distribution for ' shortfile])

vess_seg = vessels_labeled_segs;
hs_img = binary_hs_image;
skel_img = skel_img;

%%Tissue boundary 
dapi_img = dapi_img;
tissue_boundary = uint16(dapi_img>600);
tissue_thresh = tissue_boundary>0;
tissue_vol = sum(tissue_thresh(:));

%%Tissue boundary for U87-GFP tumours
% gfp_img = dapi_img;
% tissue_thresh = mat2gray(gfp_img)>0.025;
% tissue_vol = sum(tissue_thresh(:));

%%%% Vessel Diffusion Analysis
hs_img_bin = hs_img>0;
skel_img_bin = skel_img>0;

ves_segs_img = vess_seg;
vess_segs_img_bin = ves_segs_img>0;

SE_cube = strel('cube',3); 
vess_segs_bin_erode = imerode(vess_segs_img_bin,SE_cube);
vess_surface = (vess_segs_img_bin-vess_segs_bin_erode)>0;

num_hs = sum(hs_img_bin(:));
sum_vesslength = sum(skel_img_bin(:));
even_hs_dist = sum_vesslength/num_hs;
voxels_skel = regionprops3(vess_surface,'VoxelList','VoxelIdxList');
even_hs_loc = cell(size(voxels_skel,1),1);
px_per_um = px_per_um;

vessel_length_um = sum_vesslength*px_per_um;

vessel_length_mm = vessel_length_um/1000;

num_hs_even_density = round(vessel_length_mm*32.4);


%%Randomly picking # HS from vessel outline list
all_skel_voxels = cell2mat(voxels_skel.VoxelList);
rdm_pts = randperm(size(all_skel_voxels,1),num_hs);
random_hs_pts = zeros(num_hs,3);
for w = 1:num_hs
   random_hs_pts(w,:)=all_skel_voxels(rdm_pts(w),:);   
end

binary_hs_image_even = zeros(size(skel_img_bin));

for u = 1:size(random_hs_pts,1)
binary_hs_image_even(random_hs_pts(u,2),random_hs_pts(u,1),random_hs_pts(u,3))=1;
end

dist_tform_hs_even = bwdist(binary_hs_image_even);
dist_tform_hs = bwdist(hs_img_bin);


img_actual_dif_crop = dist_tform_hs.*tissue_thresh;
img_even_dif_crop = dist_tform_hs_even.*tissue_thresh;

actual_dif_img = img_actual_dif_crop;
even_dif_img = img_even_dif_crop;

%%Distances
convert_px_to_um = 100/px_per_um;

bin_vess = vess_seg>0;
dist_vess = bwdist(bin_vess);
dist_vess_less100 = dist_vess<convert_px_to_um;

tissue_boundary = tissue_thresh;
tissue_boundary = uint16(tissue_boundary)-uint16(bin_vess);

dist_vess_crop = uint16(dist_vess).*tissue_boundary;

even_vessel100_img_bin = uint16(dist_vess_less100).*tissue_boundary;
hs_vessel100_img_bin = uint16(actual_dif_img<convert_px_to_um).*tissue_boundary;
even_hs_vessel100_img_bin = uint16(even_dif_img<convert_px_to_um).*tissue_boundary;

total_tissue_vol = sum(tissue_boundary(:));
total_100fromhs = sum(hs_vessel100_img_bin(:));
total_100fromvessel = sum(even_vessel100_img_bin(:));
total_100fromhs_rando = sum(even_hs_vessel100_img_bin(:));


pcent_t100vess = total_100fromvessel/total_tissue_vol*100;
pcent_t100hs = total_100fromhs/total_tissue_vol*100;
pcent_t100hs_rando = total_100fromhs_rando/total_tissue_vol*100;

results_table = table('Size', [1 3],'VariableTypes',{'double' 'double' 'double'}, 'VariableNames',{'pcent_100um_from_all_ves_surfaces' 'pcent_100um_from_actual_ntecs' 'pcent_100um_from_random_ntecs'}); 
results_table.pcent_t100vess = pcent_t100vess;
results_table.pcent_t100hs = pcent_t100hs;
results_table.pcent_t100hs_rando = pcent_t100hs_rando;


even_along_vess = dist_vess_crop;
random_ntecs = img_even_dif_crop;
actual_ntecs = img_actual_dif_crop;

cd(save_dir)

table_name = strcat(shortfile,'_Results_Diff_voloftumour.csv');
writetable(results_table,table_name)

img_dif_ves_name = strcat(shortfile,'_evendiffalongvess.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(dist_vess_crop), img_dif_ves_name, options);  

img_actualhs_name = strcat(shortfile,'_actual_ntecs_dif.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(img_actual_dif_crop), img_actualhs_name, options);
            
img_evenhs_name = strcat(shortfile,'_random_ntecs_dif.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(img_even_dif_crop), img_evenhs_name, options);             
            


toc
end
