function [img_skel, branch_pts_image, labeled_vess_segs, vess_radius, binary_hs_image, hotspots_pts_image_dil, results_table_um_new_table, hs_diameter_table] = vessel_ntec_analysis_func(vessels_seg, dapi_img, particles, px_per_um, save_dir,sample_name)
tic
shortfile = sample_name;
display(['Analyzing vessels for ' shortfile])


%%%% Vessel post-processing, skeletonization and vessel chopping into
%%%% individual vessel segments

ves_thresh_img = vessels_seg;


num_slices = size(ves_thresh_img);
ves_thresh_bin = ves_thresh_img==1;
ves_thresh_bin = ves_thresh_bin>0;
se_opening = strel('sphere',3);
se_closing = strel('sphere',3);
ves_thresh_close = imclose(ves_thresh_bin,se_closing);
ves_thresh_erode = imerode(ves_thresh_close,se_opening);
ves_thresh_open = imdilate(ves_thresh_erode,se_opening);

ves_thresh = ves_thresh_open;

slow_img_skel = bwskel(ves_thresh,'MinBranchLength',15);
slow_img_skel = slow_img_skel>0;

img_skel = slow_img_skel;

%slow_img_skel = bwskel(ves_thresh,'MinBranchLength',15);
branch_pts = bwmorph3(slow_img_skel,'branchpoints'); 
vessel_seg_skel = slow_img_skel-branch_pts;
end_pts = bwmorph3(vessel_seg_skel,'endpoints');
label_vessel_seg_skel = bwlabeln(vessel_seg_skel);

stats = regionprops3(label_vessel_seg_skel,'Volume','VoxelIdxList');

labeled_vess_segs = im2mat(dip_growregions(uint32(label_vessel_seg_skel),[],ves_thresh,3,10,'low_first'));

labeled_vess_seg_ext_regions = im2mat(dip_growregions(uint32(labeled_vess_segs),[],[],3,2,'low_first'));

SE = strel('sphere',6);           
branch_pts_image = imdilate(branch_pts, SE);

%%Hotspot assignment
particles = particles;

dist_tform_vessel = bwdist(ves_thresh);

dist_threshold = dist_tform_vessel<5;
% SE_erode = strel('sphere',1);           
% erode_vessel = imerode(ves_thresh,SE_erode);
% 
% peri_ves_regions = (dist_threshold-erode_vessel)>0;
np_near_vess = particles.*uint16(dist_threshold);

vess_NP_Pixels = np_near_vess(np_near_vess > 0);
threshold_NP_near_ves = graythresh(vess_NP_Pixels);
treshold_np = imbinarize(np_near_vess,threshold_NP_near_ves*3);
label_hotspots = bwlabeln(treshold_np);

hotspots_analysis = regionprops3(label_hotspots,np_near_vess,'WeightedCentroid','EquivDiameter');
binary_hs_image = zeros(size(np_near_vess));

rounded_pts = round(hotspots_analysis.WeightedCentroid);

for r = 1:size(rounded_pts,1)
binary_hs_image(rounded_pts(r,2),rounded_pts(r,1),rounded_pts(r,3))=1;
end

label_hs = bwlabeln(binary_hs_image);
SE_hs = strel('sphere',5);           
hotspots_pts_image_dil = imdilate(label_hs, SE_hs);

hotspot_ves_segs = regionprops3(labeled_vess_seg_ext_regions,binary_hs_image,'VoxelValues');
array_hotspot_ves_segs = table2array(hotspot_ves_segs);

hs_per_ves_seg = zeros(size(array_hotspot_ves_segs,1),1);

for t=1:length(array_hotspot_ves_segs)
hs_per_ves_seg(t)=sum(array_hotspot_ves_segs{t},1);
end

hs_per_ves_seg_table = array2table(hs_per_ves_seg);

hs_per_ves_seg_table.Properties.VariableNames = {'Num_hotsposts_per_ves_seg'};

%%Assign pixel dimensions in um
px_per_um = px_per_um;

%%hotspot diameter
HS_diameter_table = hotspots_analysis.EquivDiameter;

results_table_dia_um_new = [HS_diameter_table(:,1).*px_per_um];
results_table_dia_um_new = array2table(results_table_dia_um_new);

results_table_dia_um_new.Properties.VariableNames = {'Hotspot_diameter_um'};
hs_diameter_table = results_table_dia_um_new;

%%Vess diameter distance transform
vess_neg = ves_thresh~=1;
vess_neg_dia = bwdist(vess_neg);

vess_radius = vess_neg_dia; 

%%%%Vessel Length using skeleton where length = volume
ves_length = regionprops3(label_vessel_seg_skel,'Volume');

ves_length.Properties.VariableNames = {'Vessel_length_px'};

%%%%Vessel Diameter from average radius mulitplied by 2
ves_rad = regionprops3(label_vessel_seg_skel,vess_neg_dia,'MeanIntensity');
ves_dia_array = (table2array(ves_rad))*2;
ves_dia_table = array2table(ves_dia_array);

ves_dia_table.Properties.VariableNames = {'Vessel_diameter_px'};

%%%%Vessel Surface area and Volume of each vessel segment
ves_sa_vol = regionprops3(labeled_vess_segs,'Volume','SurfaceArea','Centroid');
vess_sa_vol_comb = [ves_sa_vol.Volume ves_sa_vol.SurfaceArea];
vess_sa_vol_table = array2table(vess_sa_vol_comb);

vess_sa_vol_table.Properties.VariableNames = {'Vessel_Volume_px' 'Vessel_Surface_Area_px'};

%%%%Vessel Surface area:Volume
ves_sa_to_vol_ratio = ves_sa_vol.SurfaceArea./ves_sa_vol.Volume;
ves_sa_to_vol_ratio_table = array2table(ves_sa_to_vol_ratio);

ves_sa_to_vol_ratio_table.Properties.VariableNames = {'Vessel_SA_to_Vol_Ratio'};

%%%%Vessel Distance to nearest junction point 
branch_pts_Centroid = regionprops3(branch_pts,'Centroid');
dist_allpopints = pdist2(ves_sa_vol.Centroid,branch_pts_Centroid.Centroid);
ves_dist_to_branch_pt = min(dist_allpopints,[],2);
ves_dist_to_branch_pt_table = array2table(ves_dist_to_branch_pt);

ves_dist_to_branch_pt_table.Properties.VariableNames = {'Vessel_dist_to_branch_pt_px'};

%%%%Vessel Distance to nearest vessel segment, nearest 5 vessels, nearest
%%%%10 vessels 
dist_all_ves = pdist(ves_sa_vol.Centroid);
dist_all_ves_sq = squareform(dist_all_ves);

dist_all_ves_min = mink(dist_all_ves_sq,2,2);
dist_all_ves_min_5 = mink(dist_all_ves_sq,6,2);
dist_all_ves_min_10 = mink(dist_all_ves_sq,11,2);

dist_all_ves_min_5 = dist_all_ves_min_5(:,2:6);
dist_all_ves_min_10 = dist_all_ves_min_10(:,2:11);
dist_all_ves_min = dist_all_ves_min(:,2);

dist_all_ves_min_5_mean = mean(dist_all_ves_min_5,2);
dist_all_ves_min_10_mean = mean(dist_all_ves_min_10,2);

dist_all_ves_closest_table = array2table(dist_all_ves_min);
dist_all_ves_min_5_mean_table = array2table(dist_all_ves_min_5_mean);
dist_all_ves_min_10_mean_table = array2table(dist_all_ves_min_10_mean);

dist_all_ves_closest_table.Properties.VariableNames = {'Vessel_dist_to_nearest_ves_px'};
dist_all_ves_min_5_mean_table.Properties.VariableNames = {'Vessel_dist_to_nearest_5_ves_px'};
dist_all_ves_min_10_mean_table.Properties.VariableNames = {'Vessel_dist_to_nearest_10_ves_px'};

%%%%Vessel Tortuosity 
num_vessels = max(label_vessel_seg_skel(:));

dist_endpoints = zeros(num_vessels,1);

for x = 1:num_vessels
   seg_mask = label_vessel_seg_skel;
   seg_mask(seg_mask~=x) = 0;
   isolate_seg = (seg_mask>0);
   isolate_endpoints =  isolate_seg.*end_pts;
   if sum(isolate_endpoints(:)) == 0
      end_point_dist = 0; 
      dist_endpoints(x,1) = end_point_dist; 
   else 
   isolate_endpoints_label = bwlabeln(isolate_endpoints);
   end_point_coords = regionprops3(isolate_endpoints_label,'Centroid');
   end_point_dist = pdist(end_point_coords.Centroid);
   size_end_pt = size(end_point_coords.Centroid);
   if size_end_pt(1)==1
      end_point_dist = 1; 
      dist_endpoints(x,1) = end_point_dist;    
   else
   dist_endpoints(x,1) = end_point_dist;      
   end
   end
end

ves_length_array = table2array(ves_length);
tortuosity = ves_length_array./dist_endpoints;
tortuosity(tortuosity == Inf) = NaN;

tortuosity_table = array2table(tortuosity);


tortuosity_table.Properties.VariableNames = {'Vessel_tortuosity'};




%%%%Mean NP Intensity Around Vessel Segments
dapi_img = dapi_img;

bin_vess = ves_thresh>0;
dist_vess = bwdist(bin_vess);
dist_vess_less50 = dist_vess<100;


tissue_boundary = uint16(dapi_img>2500);
tissue_boundary_log = tissue_boundary>0;

dist_vess_less50_crop = dist_vess_less50.*tissue_boundary_log;
dist_vess_less50_crop = dist_vess_less50_crop>0;


tissue_boundary_vess = tissue_boundary-uint16(bin_vess);

labeled_vess_seg_ext_regions_2 = im2mat(dip_growregions(uint32(labeled_vess_segs),[],dist_vess_less50_crop,3,100,'low_first'));

np_img_crop = uint16(particles).*tissue_boundary_vess;

mean_np_int_per_seg = regionprops3(labeled_vess_seg_ext_regions_2,np_img_crop, 'MeanIntensity');
mean_np_int_per_seg.Properties.VariableNames = {'Mean_Np_Int_100um'};


%%%%Create final table

results_table = [ves_length ves_dia_table vess_sa_vol_table ves_sa_to_vol_ratio_table ves_dist_to_branch_pt_table dist_all_ves_closest_table dist_all_ves_min_5_mean_table dist_all_ves_min_10_mean_table tortuosity_table mean_np_int_per_seg hs_per_ves_seg_table];

results_table_um = table2array(results_table);

results_table_um_new = [results_table_um(:,1).*px_per_um results_table_um(:,2).*px_per_um results_table_um(:,3).*(px_per_um*px_per_um*px_per_um) results_table_um(:,4).*(px_per_um*px_per_um) results_table_um(:,5) results_table_um(:,6).*px_per_um results_table_um(:,7).*px_per_um results_table_um(:,8).*px_per_um results_table_um(:,9).*px_per_um results_table_um(:,10) results_table_um(:,11) results_table_um(:,12) ];

results_table_um_new_table = array2table(results_table_um_new);


results_table_um_new_table.Properties.VariableNames = {'Vessel_length_um' 'Vessel_diameter_um' 'Vessel_vol_um' 'Vessel_SA_um' 'Vessel_SA_to_Vol_ratio' 'Vessel_dist_to_branch_pt_um' 'Vessel_dist_to_nearest_vessel_um' 'Vessel_dist_to_nearest_5_vessel_um' 'Vessel_dist_to_nearest_10_vessel_um' 'Vessel_tortuosity' 'Mean_Np_Int_100um' 'Num_hotspots_per_ves_seg'};


%%Save Tables in save_dir

cd(save_dir)
table_name = strcat(shortfile,'_Results_Vessel_analysis.xlsx');
table_name_um = strcat(shortfile,'_Results_Vessel_analysis_um.xlsx');

writetable(results_table,table_name);
writetable(results_table_um_new_table,table_name_um);

table_name_dia = strcat(shortfile,'_Results_hotspot_dia_analysis.xlsx');

writetable(hs_diameter_table,table_name_dia);

% table_name_diff = strcat(shortfile,'_Diffusion_data.csv');
% writetable(result_test_table_diff,table_name_diff);

%%%% Write Image files to save_dir

cd(save_dir)

img_skel_name = strcat(shortfile,'_skeleton.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(slow_img_skel), img_skel_name, options);
            
branch_pts_image_name = strcat(shortfile,'_branch_pts.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(branch_pts_image), branch_pts_image_name, options);
            
labeled_vess_segs_name = strcat(shortfile,'_labeled_vess_segs.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(labeled_vess_segs), labeled_vess_segs_name, options);            

vess_radius_segs_name = strcat(shortfile,'_vess_radius.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(vess_neg_dia), vess_radius_segs_name, options);   

vess_hs_segs_name = strcat(shortfile,'_bin_hs_img.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(binary_hs_image), vess_hs_segs_name, options);  

vess_dil_hs_segs_name = strcat(shortfile,'_dilated_hs_img.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(hotspots_pts_image_dil), vess_dil_hs_segs_name, options);   
            
                  
toc



end

