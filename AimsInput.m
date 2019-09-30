% Naming conventions
% https://visualgit.readthedocs.io/en/latest/pages/naming_convention.html
% Make a class for the data or from PredictGCDMS send the data to this
% class. This class is specifically for Predict GCDMS. Can rename all the
% data here. Constructor can rename the classes. This will help reduce the
% number of variables in workspace and make it clear what variable is used
% for what. Possibly make another class that inherits from this class to.
% PredictGCDMS will create an object of this class. and will make it easier
% to read what is what.  

% Access set both setAccess and getAccess
% SetAccess relates only to changing the variable with . operator outside
% of class
% Get Access relates only to scope of getting the variables with . operator
% outside of class
classdef AimsInput
    % Create instance variable and makes them private so that other
    % function cannot access and modify them. Only memeber functions can
    % access.
    % Later change to Access. SetAccess here so I can view what object has
    % in other methods
    properties (SetAccess = private)
        sample_names = {};
        sample_labels = {};
        compensation_voltage = {};
        retention_time = {};
        intensity = {};
        watershed_label = {};
        sel_sample_index = 1; 
        bw_intensity = {};
        watershed_index = {};
        watershed_crvi = {};
        peak_table = {};
    end
    % private methods are methods that can only be used in the class. When
    % manipulating the variables it's best to have private methods
    methods 
        % Constructor must return object of class. Can change name of obj
        % parameter - (cell,cell,cell,cell)
        function obj = AimsInput(sample_names,...
                                compensation_voltage,...
                                retention_time,...
                                intensity)
            obj.sample_names = sample_names;
            obj.compensation_voltage = compensation_voltage;
            obj.retention_time = retention_time;
            obj.intensity = intensity;
            %obj.watershed_label = cell(numel(sample_names),1);
            obj.bw_intensity = cell(numel(sample_names),1);
        end
        % Mutator Methods used to recompute any imaging
        function obj = compute_all_watershed(obj,level)
            % need to compute watershed for all samples
            % Start moving in code from the script and then use private
            % methods to chunk out big concepts in the script
            % This iterates through all samples
            num_sample = numel(obj.sample_names);
            peak_tables = cell(num_sample,1);
            for i=1:num_sample
                % celldata only used for noise
                obj = obj.compute_one_watershed(i,level);
                %{
                pk_table = obj.compute_one_watershed(i);
                peak_tables{i,1} = pk_table;
                %}
            end
        end
        % This function is used to compute all peak statistics for all
        % samples.
        % Uses watershed_index to compute all statistics
        function obj = compute_peak_statistics(obj)
            obj.watershed_crvi = {};
            for sample_num=1:(size(obj.sample_names,1))
                % Go through every peak and compute peak_vol
                % temp_intensity = obj.watershed_intensity{i};
                num_of_peaks = size(obj.watershed_index{sample_num},1);
                for peak_num=1:(num_of_peaks)
                    index_list = cell2mat(obj.watershed_index{sample_num}(peak_num,1));
                    [max_cv,max_rt,max_intensity] = obj.compute_max_cv_rt_intensity(sample_num,index_list); 
                    obj.watershed_crvi{sample_num,1}(peak_num,1) = max_cv;
                    obj.watershed_crvi{sample_num,1}(peak_num,2) = max_rt;
                    obj.watershed_crvi{sample_num,1}(peak_num,3) = obj.compute_volume(sample_num,index_list);
                    obj.watershed_crvi{sample_num,1}(peak_num,4) = max_intensity;
                end
            end
        end
        % This function is used to filter out peaks that are within a
        % certain standard deviation
        function obj = filter_std_peak(obj,lower_cv,upper_cv,lower_rt,upper_rt,sd)
            % Compute noise value
            noise_index = obj.compute_index_list(lower_cv,upper_cv,lower_rt,upper_rt);
            noise = obj.intensity{obj.sel_sample_index}(noise_index);
            noise(noise<0) == 0;
            noise_sd = std(noise);
            if (sd > 0)
                noise_sd = noise_sd * sd;
            end

            % Go through every sample and apply standard deviation to all
            % peaks
            for i=1:(size(obj.sample_names,1))
                % volume is in third column
                peak_intensities = obj.watershed_crvi{i}(:,4); 
                peak_num = (1:1:size(peak_intensities,1))';
                peak_num_intensities = [peak_num,peak_intensities];
                filter_peak_num = obj.ret_noise_peak(peak_num_intensities,noise_sd)';
                if (size(filter_peak_num,1))
                    [watershed_label,watershed_index] = obj.filter_watershed_label_index(obj.watershed_label{i},obj.watershed_index{i},filter_peak_num);
                    obj.watershed_label{i} = watershed_label;
                    obj.watershed_index{i} = watershed_index;
                end
            end   
        end
        function obj = generate_peak_table(obj,s1_index,s2_index,cv_tolerance,rt_tolerance)
            % Split samples into healthy and unhealthy
            % Then choose a sample that will be representative of the
            % sample
            % For every sample store the cv and rt ranges
            % For reference choose the one with the most amount of peaks
            %s1_index = 1;
            %s2_index = 2;
            %cv_tolerance = 1; % 1
            %rt_tolerance = 10; % 10
            sample_labels = obj.sample_labels;
            watershed_crvi = obj.watershed_crvi;
            unique_labels = string(unique(sample_labels));
            
            s1 = watershed_crvi{s1_index}(:,1:2);
            peak_num = (1:1:size(s1,1))';
            s1 = [peak_num,s1];
            
            s2 = watershed_crvi{s2_index}(:,1:2);
            peak_num = (1:1:size(s2,1))';
            s2 = [peak_num,s2];
            
            s1_s2_overlap = [];
            for i=1:size(s1,1)
                curr_s1 = s1(i,:);
                s2_overlap = obj.ret_overlap(curr_s1,s2,cv_tolerance,rt_tolerance);
                if (size(s2_overlap,1) == 0)
                    continue;
                end
                % if there are multiple peaks that overlap then we choose
                % the one that has the closest euclidean distance
                if(size(s2_overlap,1)>1)
                    % Compute euclidean distance across all overlap
                    euclidean_dist = ((curr_s1(1,2)-s2(s2_overlap,2)).^2 + (curr_s1(1,3)-s2(s2_overlap,3)).^2).^(1/2);
                    [min_val,min_index] = min(euclidean_dist);
                    s2_overlap = s2_overlap(min_index);
                end
                add_overlap = [i,s2_overlap];
                s1_s2_overlap = [s1_s2_overlap;add_overlap];
            end
            
            % use s1_s2_overlap and if any duplicate in the second column
            % then remove them from s1_s2_overlap. Those will form their
            % own column. 
            % For the non-overlap just use the peak cv and rt as the 
            % For over-lap take the average of the peak cv and rt. should
            % be used.
            [~,unique_index,~] = unique(s1_s2_overlap(:,2));
            dup_index = setdiff((1:size(s1_s2_overlap, 1))', unique_index);
            dup_val = s1_s2_overlap(dup_index,:);
            % may have multiple duplicate values
            dup_val = unique(dup_val(:,2));
            non_overlap_filter = ~ismember(s1_s2_overlap(:,2),dup_val);
            s1_s2_overlap = s1_s2_overlap(non_overlap_filter,:);
            
            % Construct the peak table. The columns needs to store the max
            % cv/rt coordinates. For the overlap take th average between
            % both peaks.
            peak_table_cv_rt = [];
            for i=1:size(s1_s2_overlap,1)
                cv_rt = [s1(s1_s2_overlap(i,1),2:3);s2(s1_s2_overlap(i,2),2:3)];
                peak_avg_cv_rt = mean(cv_rt);
                peak_table_cv_rt = [peak_table_cv_rt;peak_avg_cv_rt];
            end
            
            s1_non_overlap_index = ~ismember(s1(:,1),s1_s2_overlap(:,1));
            s1_cv_rt = s1(s1_non_overlap_index,2:3);
            peak_table_cv_rt = [peak_table_cv_rt;s1_cv_rt];
            
            s2_non_overlap_index = ~ismember(s2(:,1),s1_s2_overlap(:,2));
            s2_cv_rt = s2(s2_non_overlap_index,2:3);
            peak_table_cv_rt = [peak_table_cv_rt;s2_cv_rt];
            
            % go through all samples and see is the peak exists. If it does
            % then log it into peak table
            num_samples = size(obj.sample_names,1);
            num_peaks = size(peak_table_cv_rt,1);
            peak_table = zeros(num_samples,num_peaks);
            
            for sample = 1:num_samples
                curr_peak_crvi = obj.watershed_crvi{sample};
                curr_peak_crv = curr_peak_crvi(:,1:3);
                low_up_cv = [(peak_table_cv_rt(:,1) - cv_tolerance),(peak_table_cv_rt(:,1) + cv_tolerance)];
                low_up_rt = [(peak_table_cv_rt(:,2) - rt_tolerance),(peak_table_cv_rt(:,2) + rt_tolerance)];
                for peak = 1:num_peaks
                    temp_peak_crv = curr_peak_crv(curr_peak_crv(:,1) >= low_up_cv(peak,1),:);
                    temp_peak_crv = temp_peak_crv(temp_peak_crv(:,1) <= low_up_cv(peak,2),:);
                    temp_peak_crv = temp_peak_crv(temp_peak_crv(:,2) >= low_up_rt(peak,1),:);
                    temp_peak_crv = temp_peak_crv(temp_peak_crv(:,2) <= low_up_rt(peak,2),:);
                    num_peak_sat = size(temp_peak_crv,1);
                    if (num_peak_sat == 0)
                        continue;   
                    end
                    if (num_peak_sat > 1)
                        % filter out peaks and choose the one that is
                        % closest
                        euclidean_dist = ((peak_table_cv_rt(peak,1)-temp_peak_crv(:,1)).^2 + (peak_table_cv_rt(peak,2)-temp_peak_crv(:,2)).^2).^(1/2);
                        [min_val,min_index] = min(euclidean_dist);
                        temp_peak_crv = temp_peak_crv(min_index,:);
                    end
                    peak_table(sample,peak) = temp_peak_crv(1,3);                  
                end
            end
            obj.peak_table = peak_table;
        end
        
        
        
        % Access Methods to retrieve instance variables the type of data return is
        % on leftside of equal sign
        function obj = set_bw_intensity(obj,index,bw_intensity)
            %bw_intensity = imcomplement(bw_intensity);
            obj.bw_intensity{index} = int8(bw_intensity);
        end 
        function obj = set_sel_sample_index(obj,index)
            obj.sel_sample_index = index;
        end
        function obj = set_sample_labels(obj,new_labels)
            obj.sample_labels = new_labels;
        end
        function cell = get_sample_names(obj)
            cell = obj.sample_names;
        end
        function mat_double = get_cv(obj,index)
            mat_double = obj.compensation_voltage{index};
        end
        function mat_double = get_rt(obj,index)
            mat_double = obj.retention_time{index};
        end
        function mat_double = get_intensity(obj,index)
            mat_double = obj.intensity{index};
        end
        function mat_double = get_watershed_intensity(obj,index)
            mat_double = obj.watershed_label{index};
        end
        function int = get_sel_sample_index(obj)
            int = obj.sel_sample_index;
        end
        function mat_double = get_bw_intensity(obj,index)
            mat_double = obj.bw_intensity{index};
        end
        function mat_double = get_watershed_label(obj,index)
            mat_double = obj.watershed_label{index};
        end
        function mat_double = get_peak_table(obj)
            mat_double = obj.peak_table;
        end
        
    end
    % private mutator methods
    methods (Access = private)
        % computes watershed on one sample 
        % parameters - (double,double,double)
        function obj = compute_one_watershed(obj,index,level)
            cv = obj.compensation_voltage{index};
            rt = obj.retention_time{index};
            intensity = obj.intensity{index};
            disk_radius = 12;
            
            
            calc_intensity = intensity;
            % Sets all negative inensity to zero
            calc_intensity(calc_intensity <= 0) = 0;
            % Turn intensity to grayscale image
            max_intensity = max(max(calc_intensity));
            calc_intensity = mat2gray(calc_intensity,[0,max_intensity]);
            calc_intensity = obj.top_hat_filter(calc_intensity,disk_radius);
            [bin_calc_intensity, bw_intensity] = obj.binarize(calc_intensity,level);
            % save binarized data for graphing
            obj = obj.set_bw_intensity(index,bw_intensity);
            
            watershed_label = obj.alg_watershed(bin_calc_intensity);
            [watershed_label,watershed_index] = obj.filter_background(index,watershed_label);
            
            % save watershed intensity
            obj.watershed_label{index,1} = watershed_label;
            obj.watershed_index{index,1} = watershed_index;
        end
        % Other shapes of structuring element
        % https://www.mathworks.com/help/images/ref/strel.html
        function mat_double_filtered = top_hat_filter(obj,intensity,disk_radius)
            se = strel('disk',disk_radius); 
            intensity_filtered = imtophat(intensity,se);
            mat_double_filtered = intensity_filtered;
        end
        function [mat_double_binarized,bw_intensity] = binarize(obj,intensity,levels)
            % improve the constrast of the image to help show the peaks
            % more
            adj_intensity = imadjust(intensity);
            bin_intensity = im2bw(adj_intensity, levels);
            bw_intensity = bin_intensity;
            % Make bakcground -inf and peaks 1
            log_background = ~bin_intensity;
            bin_intensity = -bwdist(log_background);
            bin_intensity(log_background) = -Inf;
            mat_double_binarized = bin_intensity; 
        end
        function mat_double_watershed = alg_watershed(obj,intensity)
            mat_double_watershed = watershed(intensity);
        end
        % Watershed algorithm will sometimes make the background one giant
        % peak. Sometimes it will separate the background into separate
        % giant peaks
        % In order to fix this, take the black and white intensity image
        % before it goes into watershed.
        % Take the balck and white intensity image background and overlap
        % the backgroundo onto watershed image to bring back the background
        % to watershed algoritm
        function [watershed_label,watershed_index] = filter_background(obj,index,watershed_label)
            % should use bwconnocomp to find the super large peaks that
            % cannot be possible. Can also filter out peaks that are small
            % background is 0 in bw_intensity so flip the 1 and 0s
            % filter according to size
            max_peak_size = 100;
            cc = bwconncomp(watershed_label);
            % Gets index of all peaks
            watershed_index = cc.PixelIdxList';
            % Get size of all pixel size of all peaks
            peak_size = cellfun('size',watershed_index,1);
            % filter out peaks above certain range
            filter_peak_num = find(peak_size > max_peak_size);
            [watershed_label,watershed_index] = obj.filter_watershed_label_index(watershed_label,watershed_index,filter_peak_num);
            
            % This code should be deleted but may be useful later
            % Apply filter to watershed_label 
            % bool_filter = ismember(watershed_label,filter_peak_num);
            % watershed_label(bool_filter) = 0;
            % once you filter out max peaks need to reassign numbers to
            % peaks because map doesn't match
            
        end
        function [watershed_label,watershed_index] = filter_watershed_label_index(~,watershed_label,watershed_index,filter_peak_num)
            watershed_index(filter_peak_num) = [];
            % once peak indexes are removed must repaint watershed and put
            % the intensity labels back in
            watershed_label = uint8(zeros(size(watershed_label)));
            for i = 1:size(watershed_index,1)
                watershed_label(watershed_index{i}) = i;
            end
        end
        function mat_uint_list = compute_index_list(obj,l_cv,u_cv,l_rt,u_rt)
            % Pull out the sample we need to get noise from
            index = obj.sel_sample_index;
            cv = obj.compensation_voltage{index};
            rt = obj.retention_time{index};
            
            bool_l_cv = cv<=l_cv;
            bool_u_cv = cv>=u_cv;
            bool_l_rt = rt<=l_rt;
            bool_u_rt = rt>=u_rt;
            l_cv_index = max(find(bool_l_cv));
            u_cv_index = min(find(bool_u_cv));
            l_rt_index = max(find(bool_l_rt));
            u_rt_index = min(find(bool_u_rt));
            
            bool_intensity = obj.intensity{index};
            bool_intensity = zeros(size(bool_intensity));
            bool_intensity(l_rt_index:u_rt_index,l_cv_index:u_cv_index) = 1;
            mat_uint_list = uint32(find(bool_intensity));
        end
        function double_peak_vol = compute_volume(obj,index,index_list)
            
            cv = obj.compensation_voltage{index};
            rt = obj.retention_time{index};
            temp_intensity = obj.intensity{index};
            
            total_len = abs(abs(cv(2))-abs(cv(1)));
            total_width = abs(abs(rt(2))-abs(rt(1)));
            %{
            % for debugging
            for i=1:size(index_list,1)
                if (index_list(i,1) >= 305000)
                    disp('a');
                end
            end
            %}
            heights = temp_intensity(index_list);
            total_height = sum(sum(heights));
            double_peak_vol = total_len * total_height * total_width;
        end
        function [double_max_cv,double_max_rt,double_max_intensity] = compute_max_cv_rt_intensity(obj,sample_num,index_list)
            % if an error occurs in this function it could be due to too
            % multiple indexes with max intensity. Code in something that
            % will randomly choose a max index
            temp_intensity = obj.intensity{sample_num};
            [max_intensity,I] = max(temp_intensity(index_list));
            max_index = index_list(I);
            num_rows = size(temp_intensity,1);
            num_cols = size(temp_intensity,2);
            col = ceil(max_index./num_rows);
            row = rem(max_index,num_rows);
            % covers the edges cases
            if ((col == 0) & (row == 0))
                col = 1;
                row = 1;
            end
            if (row ==0)
                row = num_cols;
            end
            
            double_max_rt = obj.retention_time{sample_num}(row,1);
            double_max_cv = obj.compensation_voltage{sample_num}(col,1);
            double_max_intensity = max_intensity;
        end

        function list = ret_noise_peak(obj,peak_num_intensities,noise_sd)
            % base case
            if(size(peak_num_intensities,1) == 1)
                if(peak_num_intensities(1,2) < noise_sd)
                    list = peak_num_intensities(1,1);
                else
                    list = [];
                end
                return;
            end
            % recursive case
            list = [obj.ret_noise_peak(peak_num_intensities(1:end-1,:),noise_sd),obj.ret_noise_peak(peak_num_intensities(end,:),noise_sd)];
            return;
        end
        
        function list = ret_overlap(obj,s1,s2_list,cv_tol,rt_tol)
            % Base case
            if (size(s2_list,1) == 1);
                lower_cv = s1(1,2)-cv_tol;
                upper_cv = s1(1,2)+cv_tol;
                lower_rt = s1(1,3)-rt_tol;
                upper_rt = s1(1,3)+rt_tol;
                cv = s2_list(1,2);
                rt = s2_list(1,3);
                if ((cv > lower_cv && cv < upper_cv) && (rt > lower_rt && rt < upper_rt))
                    list = s2_list(1,1);
                else
                    list = [];
                end
                return;
            end
            % Recursive list
            list = [obj.ret_overlap(s1,s2_list(1:end-1,:),cv_tol,rt_tol);obj.ret_overlap(s1,s2_list(end,:),cv_tol,rt_tol)];
            return;
        end

    end
end