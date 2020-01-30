clear all;
clc; close all
dbstop if error
%----------------------------------------------------------------------------------------
% Structure Parsing
% Lower stem (LS), bulge (B), upper stem (US), Loop (L), no-hairpin (N)
%----------------------------------------------------------------------------------------
global RBP  RBPbs iRBP ifig
RBP = {'PP7' 'MS2' 'Qb'}; RBPbs = {'PP7' 'MS2' 'Qb' 'Random'}; ifig = 1;
WTseq = {'uaaggaguuuauauggaaacccuua','acaugaggauuacccaugu','augcaugucuaagacagcau'};
%     IS_workComputer =0;

% if IS_workComputer
%     Path{1} = 'C:\Users\katznoa\Documents\Dropbox_temp\';
%     Path{2} = 'C:\Users\katznoa\Documents\Dropbox_temp\';
% else
%     Path{1} = 'C:\Users\Noa\Dropbox\matlab\NGS\';
%     Path{2} = 'C:\Users\Noa\Dropbox\Lab\Binding Assay\NGS\Results\Yaron\files sent\';
% end

% load([Path{1} 'BS_linear_final_v2.mat']);

% cd(Path{2});
% filename = 
file = fopen('parFile.txt','r');
AA = textscan(file,'%s','Delimiter' ,'\n');
bsStruct = {AA{1,1}(2:2:length(AA{1,1}))};
bsSeq =  {AA{1,1}(1:2:length(AA{1,1}))};
L = length(bsSeq{1,1}{1,1});
for i = 1:length(bsStruct{1,1})
   temp = bsStruct{1,1}{i,1}; 
   bsStruct{1,1}{i,1} = temp(1:L);
end

BUG= 0;
idxs=[0]; %the initialize with 0 is for the python be able to read the csv even if there are no bugs
for iRBP = 1
  
    % extract loop sequences by the string '(' '.'xN ')'
    match_loop{iRBP} = regexp(lower(bsStruct{iRBP}),'\([.]{2,}+\)'); % find first (, than dots at least two times, than )
    match_loop_len{iRBP} = regexp(lower(bsStruct{iRBP}),'\([.]{2,}+\)','match'); % find first (, than dots at least two times, than )

    % extract forward stem sequences, the first part of hairpin (((
    match_stem_f{iRBP} = regexp(lower(bsStruct{iRBP}),'\(+[.]{1,1}');
    % find a number of (, than only one dot.
    % first indice- lower stem. second indice- upper stem. if only one-
    % there's no bulge.

    % extract reverse stem sequences, the second part of hairpin )))
    match_stem_r{iRBP} = regexp(lower(bsStruct{iRBP}),'\)+'); %,'tokenExtents');
    % find a number of (, than only one dot.
    % first indice- lower stem. second indice- upper stem. if only one-
    % there's no bulge.
    match_bulge_f{iRBP} = regexp(lower(bsStruct{iRBP}),'\([.]+\(');
    match_bulge_r{iRBP} = regexp(lower(bsStruct{iRBP}),'\)[.]+\)');
end

% write the sequences for the different parts of the binding site
bug_count = 0;
for iRBP = 1
    for iSeq = 1:size(bsStruct{iRBP},1)
        if iSeq ==1
            prevlen=0;
        else
            prevlen = length(struct_mat{iRBP});
        end
        bs_len =L;% length(bsSeq{iRBP}{iSeq});
        struct_mat{iRBP}{iSeq} = zeros(5,bs_len);

        % if identify loop (hairpin sturctre)
        if  size(match_loop{iRBP}{iSeq}) == 1
            loopInds = match_loop{iRBP}{iSeq}+1; loop_len = size(match_loop_len{iRBP}{iSeq}{1},2)-2;
            struct_mat{iRBP}{iSeq}(4,loopInds:loopInds+loop_len-1) = ones(1,loop_len); % Loop

            % do we have two stems (upper and lower) for forward and less than 4 for reverse?
            if size(match_stem_f{iRBP}{iSeq},2) == 2 && size(match_stem_r{iRBP}{iSeq},2) < 4
                bulge_f_len = match_stem_f{iRBP}{iSeq}(2) - (match_bulge_f{iRBP}{iSeq}(1)+1);
                ind_bulge_f = match_bulge_f{iRBP}{iSeq}(1)+1;

                stemIndLs_f = [match_stem_f{iRBP}{iSeq}(1) match_stem_f{iRBP}{iSeq}(2)-1-bulge_f_len];
                stemIndUs_f = [match_stem_f{iRBP}{iSeq}(2) match_loop{iRBP}{iSeq}];
                stemL_len = stemIndLs_f(2) - stemIndLs_f(1)+1;
                stemU_len = stemIndUs_f(2) - stemIndUs_f(1)+1;
                struct_mat{iRBP}{iSeq}(2,ind_bulge_f:ind_bulge_f+bulge_f_len-1) = ones(1,bulge_f_len);

                struct_mat{iRBP}{iSeq}(1,stemIndLs_f(1):stemIndLs_f(2)) = ones(1,stemL_len); % LS
                struct_mat{iRBP}{iSeq}(3,stemIndUs_f(1):stemIndUs_f(2)) = ones(1,stemU_len); % US

                % if stem doesn't start from beginning, place 1's in N (no-binding site)
                if stemIndLs_f(1) ~= 1
                    struct_mat{iRBP}{iSeq}(5,1:stemIndLs_f(1)-1) = ones(1,stemIndLs_f(1)-1); % from base 1 until start of hairpin
                    find_dot_after_stemR = regexp(lower(bsStruct{iRBP}{iSeq}),')[.]');
                    if length(find_dot_after_stemR) > 1
                        find_dot_after_stemR = max(find_dot_after_stemR);
                    end
                    ind_noStruct_r = find_dot_after_stemR+1;
                    struct_mat{iRBP}{iSeq}(5,ind_noStruct_r:end) = ones(1,bs_len-ind_noStruct_r+1);
                end

                if isempty(stemU_len)
                    'NNN';
                end

                % no bulge in reverse strand
                if size(match_stem_r{iRBP}{iSeq},2) == 1
                    stemInd_r = match_stem_r{iRBP}{iSeq}(1);                   
                    struct_mat{iRBP}{iSeq}(3,stemInd_r:stemInd_r+stemU_len-1) = ones(1,stemU_len);
                    struct_mat{iRBP}{iSeq}(1,stemInd_r+stemU_len:stemInd_r+stemU_len+stemL_len-1) = ones(1,stemL_len);
                else

                    bulge_r_len = match_stem_r{iRBP}{iSeq}(2)  - match_bulge_r{iRBP}{iSeq}-1;
                    stemIndUs_r = match_stem_r{iRBP}{iSeq}(1); stemIndLs_r = match_stem_r{iRBP}{iSeq}(2);
                    ind_bulge_r = match_bulge_r{iRBP}{iSeq}+1;
                    stemL_len_r = stemL_len + stemU_len - length(stemIndUs_r:ind_bulge_r-1);

                    try
                        struct_mat{iRBP}{iSeq}(2,ind_bulge_r:ind_bulge_r+bulge_r_len-1) = ones(1,bulge_r_len);
                        struct_mat{iRBP}{iSeq}(3,stemIndUs_r:ind_bulge_r-1) = ones(1,length(stemIndUs_r:ind_bulge_r-1));
                        struct_mat{iRBP}{iSeq}(1,stemIndLs_r:stemIndLs_r+stemL_len_r-1) = ones(1,stemL_len_r);
                    catch
                        BUG=BUG+1;
%                             disp(bsSeq{1,1}(iSeq))
%                             disp(bsStruct{1,1}(iSeq))
%                         idxs = [idxs,iSeq];
                        continue
                    end %end try
                end
            % No bulge, LS = 0, US = stem
            else
                if size(match_stem_f{iRBP}{iSeq},2) == 1 && size(match_stem_r{iRBP}{iSeq},2) < 4
                    stemIndUs_f = [match_stem_f{iRBP}{iSeq}(1) match_loop{iRBP}{iSeq}];
                    stem_len = stemIndUs_f(2)-stemIndUs_f(1)+1;
                    % place 1's in Upper stem only
                    struct_mat{iRBP}{iSeq}(3,stemIndUs_f(1):stemIndUs_f(2)) = ones(1,stem_len); % US

                    % if stem doesn't start from beginning, place 1's in N (no-binding site)
                    if stemIndUs_f(1) ~= 1
                        struct_mat{iRBP}{iSeq}(5,1:stemIndUs_f(1)-1) = ones(1,stemIndUs_f(1)-1); % from base 1 until start of hairpin
                        find_dot_after_stemR = regexp(lower(bsStruct{iRBP}{iSeq}),')[.]');
                        if length(find_dot_after_stemR) > 1
                            find_dot_after_stemR = max(find_dot_after_stemR);
                        end
                        ind_noStruct_r = find_dot_after_stemR+1;
                        struct_mat{iRBP}{iSeq}(5,ind_noStruct_r:end) = ones(1,bs_len-ind_noStruct_r+1);
                    end

                    % no bulge in reverse strand
                    if size(match_stem_r{iRBP}{iSeq},2) == 1
                        stemInd_r = match_stem_r{iRBP}{iSeq}(1); 
                        struct_mat{iRBP}{iSeq}(3,stemInd_r:stemInd_r+stem_len-1) = ones(1,stem_len);
                    else
                        bulge_r_len = match_stem_r{iRBP}{iSeq}(2)  - match_bulge_r{iRBP}{iSeq}-1;
                        stemIndUs_r = match_stem_r{iRBP}{iSeq}(1); stemIndLs_r = match_stem_r{iRBP}{iSeq}(2);
                        ind_bulge_r = match_bulge_r{iRBP}{iSeq}+1;
                        stemU_r_len = length(stemIndUs_r:ind_bulge_r-1);
                        stemL_r_len = stem_len - stemU_r_len;
                        try
                            struct_mat{iRBP}{iSeq}(2,ind_bulge_r:ind_bulge_r+bulge_r_len-1) = ones(1,bulge_r_len);
                            struct_mat{iRBP}{iSeq}(3,stemIndUs_r:ind_bulge_r-1) = ones(1,stemU_r_len);
                            struct_mat{iRBP}{iSeq}(1,stemIndLs_r:stemIndLs_r+stemL_r_len-1) = ones(1,stemL_r_len);
                        catch
                            BUG=BUG+1;
%                             idxs = [idxs,iSeq];
                            continue
                        end
                    end

                % more than two stems
                else
                    struct_mat{iRBP}{iSeq}(5,1:end) = ones(1,bs_len); % US
                    struct_mat{iRBP}{iSeq}(1:4,1:end) = zeros(4,bs_len); % US
               end
            end
        % no loop detected, meaning no binding site structure
        else
            struct_mat{iRBP}{iSeq}(5,1:end) = ones(1,bs_len); % US
            struct_mat{iRBP}{iSeq}(1:4,1:end) = zeros(4,bs_len); % US
        end

        %         % if still some are missing, place 1's in N (no-binding site).
        %         probably due to dots after stem reversed and not forward
        if sum(sum(struct_mat{iRBP}{iSeq},1)) ~= bs_len
            ind_empty = find(sum(struct_mat{iRBP}{iSeq},1) == zeros(1,bs_len));
            struct_mat{iRBP}{iSeq}(5,ind_empty) = ones(1,length(ind_empty));
        end

        % bug check
        if ~isequal(sum(struct_mat{iRBP}{iSeq},1), ones(1,bs_len))
            struct_mat{iRBP}{iSeq}(5,1:end) = ones(1,bs_len); % US
            struct_mat{iRBP}{iSeq}(1:4,1:end) = zeros(4,bs_len); % US
            bug_count = bug_count + 1;
        end
    end
    struct_mat{iRBP} = struct_mat{iRBP}';
    if length(struct_mat{iRBP}) == prevlen
       idxs = [idxs,iSeq]; 
    end
end

save( 'structure_matrix_tmp.mat','struct_mat');
csvwrite('Bugs.csv',idxs-1)    
