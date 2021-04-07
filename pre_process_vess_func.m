function [pre_vessels] = pre_process_vess_func(vessels,save_dir,sample_name)

shortfile = sample_name;
display(['Pre-processing ' shortfile])


%Pre-processing of blood vessels channel
vessels_single = single(vessels);
vessels_local = vessels_single./(0.2*max(vessels_single(:))+(im2mat(gaussf(vessels_single,size(vessels_single,1)/10))));
vessels_loglocal = mat2im(log(vessels_local+0.1));

vessels_processed = uint16(stretch(vessels_loglocal));
pre_vessels = vessels_processed;

cd(save_dir)


img_vess_name = strcat(shortfile,'_pre_processed_vessels.tif');
clear options;
            options.overwrite = true;
            options.compress = 'lzw';
            saveastiff(uint16(pre_vessels), img_vess_name, options);

end

