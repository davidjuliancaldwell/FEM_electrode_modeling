%%  Code blocks from ReconstructCorticalSurface()
%   The functions below have been extracted from our main recon MATLAB
%   function which covers the entire recon process.
%
%   This code was written by previous graduate students (I believe a
%   combination of Tim Blakely, Kai Miller, and Miah Wander, but mostly
%   T. Blakely and K. J. Miller).
%
%   The functions below are called once the pial surface has been formed
%   through FreeSurfer and electrode coordinates have been hand selected in
%   BioImage Suite using the CT scan co-registered to MR space
%
%   I have added some comments throughout to hopefully provide context and
%   clarify some things
%
%   Nile Wilson, 2018.04.03

%%  Generate surface from MRI
%   This function generates the hollow brain surface (brainHull) using the
%   volumetric MRI brain. This also generates transformMat, which is used
%   later in ProjectElectrodes()

function GenSurfFromMRI(fullName, destName, threshold)
    brain_info = spm_vol(fullName);
    brainVol = double(brain_info.private.dat);
    fprintf('  Creating surface\n');
    transformMat = brain_info.mat;
    temp1 = isosurface(brainVol,threshold);
    temp2 = isonormals(smooth3(brainVol),temp1.vertices);

    cortex = temp1;
    cortex.vertexnormals  = temp2;
    cortex.facevertexcdata = zeros(size(cortex.vertices,1),1);
    rotMat = brain_info.mat(1:3,1:3)';

    newRot = [[0 1 0];[1 0 0];[0 0 1]];
    cortex.vertices = cortex.vertices * newRot;
    cortex.vertexnormals = cortex.vertexnormals * newRot;

    cortex.vertices=(cortex.vertices*transformMat(1:3,1:3)')+repmat(transformMat(1:3,4),1,length(cortex.vertices))'; 
    cortex.vertexnormals=(cortex.vertexnormals*rotMat);

    cortex.vertexnormals=cortex.vertexnormals ./ repmat(sqrt(sum((cortex.vertexnormals .* cortex.vertexnormals),2)),1,3);

    fprintf('  Loading volume for hull generation\n');
    brain_info=spm_vol(fullName); 
    [brainHull]=spm_read_vols(brain_info);

    %threshold like the isosurf does
    brainHull = (brainHull > threshold);

    %create the convex hull
    fprintf('  Creating the convex hull\n');
    for i=1:size(brainHull,3)
        [x y] = find(brainHull(:,:,i) > 0);
        if length(x) < 3 || length(y) < 3
            continue;
        end
        if length(unique(x)) < 2 || length(unique(y)) < 2
            continue;
        end
        out = convhull(x,y);
        brainHull(:,:,i) = roipoly(brainHull(:,:,i), y(out),x(out));
    end

    %blur the convex hull a little
    brainHull = sm_filt(brainHull,1);
    %get every pixel that is above threshold into the hull
    brainHull = (brainHull > 0);
    %hollow out the hull
    brainHull = hollow_brain(brainHull);

    transformMat = brain_info.mat;

    fprintf('  Saving results\n');
    eval(sprintf('save %s cortex brainHull transformMat;',destName));
end

%%  Generate needed surfaces
%   This function creates the object and .mat files for the surface, which
%   will be used in ProjectElectrodes()

function GenerateNeededSurfaces(baseDir)

    for side = {'lh','rh','both'}
        for type = {'hires','printable','lowres'}
            fileName = ['surf/' getenv('recon_patientCode') '_cortex_' side{:} '_' type{:}  '.mat'];
            isThere = TestFileExists(fileName,baseDir);
            
            if isThere == 1
                continue;
            end
            fprintf('Generating surface file ''%s''...\n', fileName);
            
            switch side{:}
                case 'both'
                    switch type{:}
                        case 'hires'
                            fprintf('  Combining left and right hi-res OBJ\n');
                            
                            meshLabExe = fullfile(myGetenv('matlab_devel_dir'), 'external', 'MeshLab', 'meshlabserver.exe');
                            inLeft = [baseDir 'surf\obj\lh.obj'];
                            inRight = [baseDir 'surf\obj\rh.obj'];
                            bothFile = [baseDir 'surf\obj\both.obj'];
                            scriptFile = fullfile(myGetenv('matlab_devel_dir'), 'Visualization', 'Recon', 'merge_lr.mlx');
                            [resulta result] = system(sprintf('%s -i %s %s -s %s -o %s -om vn', meshLabExe, inLeft, inRight,scriptFile, bothFile));
                            clear resulta result
                            
                            fprintf('  Loading combined hi-res OBJ (this may take a few minutes)\n');
                            [cortex.vertices cortex.faces cortex.normals] = load_tobj(bothFile);
                            cortex.facevertexcdata = zeros(size(cortex.vertices,1),1);

                            fprintf('  Saving\n');
                            eval(['save ' fileName ' cortex']);
                            clear cortex;
                        case 'lowres'
                            fprintf('  Combining left and right low-res OBJ\n');
                            
                            meshLabExe = [myGetenv('matlab_devel_dir') '\external\MeshLab\meshlabserver.exe'];
                            inLeft = [baseDir 'surf/obj/lh_lowres.obj'];
                            inRight = [baseDir 'surf/obj/rh_lowres.obj'];
                            bothFile = [baseDir 'surf/obj/both_lowres.obj'];
                            scriptFile = [myGetenv('matlab_devel_dir') 'Visualization\Recon\merge_lr.mlx'];
                            
                            [resulta result] = system(sprintf('%s -i %s %s -s %s -o %s -om vn', meshLabExe, inLeft, inRight,scriptFile, bothFile));
                            clear resulta result
                            
                            fprintf('  Loading combined hi-res OBJ (this may take a few minutes)\n');
                            [lowres_cortex.vertices lowres_cortex.faces lowres_cortex.normals] = load_tobj(bothFile);
                            lowres_cortex.facevertexcdata = zeros(size(lowres_cortex.vertices,1),1);

                            fprintf('  Saving\n');
                            eval(['save ' fileName ' lowres_cortex']);
                            clear lowres_cortex;
                        case 'printable'
                            
                            fprintf('  Exporting lh printable OBJ\n');
                            load([baseDir 'surf/' getenv('recon_patientCode') '_cortex_lh_printable.mat']);
                            write_obj([baseDir 'surf/obj/lh_printable.obj'],cortex.vertices,cortex.faces);
                            fprintf('  Exporting rh printable OBJ\n');
                            load([baseDir 'surf/' getenv('recon_patientCode') '_cortex_rh_printable.mat']);
                            write_obj([baseDir 'surf/obj/rh_printable.obj'],cortex.vertices,cortex.faces);
                            
                            fprintf('  Combining left and right printable OBJ\n');
                            
                            meshLabExe = [myGetenv('matlab_devel_dir') '\external\MeshLab\meshlabserver.exe'];
                            inLeft = [baseDir 'surf/obj/lh_printable.obj'];
                            inRight = [baseDir 'surf/obj/rh_printable.obj'];
                            bothFile = [baseDir 'surf/obj/both.obj'];
                            scriptFile = [myGetenv('matlab_devel_dir') 'Visualization\Recon\merge_lr.mlx'];
                            [resulta result] = system(sprintf('%s -i %s %s -s %s -o %s -om vn', meshLabExe, inLeft, inRight,scriptFile, bothFile));
                            clear resulta result
                            
                            fprintf('  Loading combined printable OBJ (this may take a few minutes)\n');
                            [cortex.vertices cortex.faces cortex.normals] = load_tobj(bothFile);
                            cortex.facevertexcdata = zeros(size(cortex.vertices,1),1);

                            fprintf('  Saving\n');
                            eval(['save ' fileName ' cortex']);
                            clear cortex;
                    end

                otherwise
                    switch type{:}
                        case 'hires'
                            threshold = 0.1;
                            GenSurfFromMRI([baseDir '/mri/' side{:} '.dpial.ribbon.nii'], [baseDir fileName], threshold);
                        case 'printable'
                            threshold = -1;
                            GenSurfFromMRI([baseDir '/mri/' side{:} '.dpial.ribbon.nii'], [baseDir fileName], threshold);
                        case 'lowres'
                            %Note: needs to be run AFTER hires, since it assumes
                            %the hires file exists
                            clear cortex;
                            [a b c] = mkdir([baseDir 'surf/obj']);
                            clear a b c;
                            fprintf('  Loading hires surface\n');
                            eval(sprintf('load %s;', ['surf/' getenv('recon_patientCode') '_cortex_' side{:} '_hires.mat']));

                            fprintf('  Exporting to OBJ\n');
                            write_obj([baseDir 'surf/obj/' side{:} '.obj'],cortex.vertices,cortex.faces);
                            fprintf('  Executing meshlab surface decimation\n');

                            % HARDCODED
                            meshLabExe = [myGetenv('matlab_devel_dir') '\external\MeshLab\meshlabserver.exe'];
                            inputFile = [baseDir 'surf/obj/' side{:} '.obj'];
                            scriptFile = [myGetenv('matlab_devel_dir') '\Visualization\Recon\create_lowpoly.mlx'];
                            outputFile = [baseDir 'surf/obj/' side{:} '_lowres.obj'];
                            [resulta result] = system(sprintf('%s -i %s -s %s -o %s -om vn', meshLabExe, inputFile, scriptFile, outputFile));
                            clear resulta result

                            fprintf('  Loading low poly OBJ\n');
                            [lowres_cortex.vertices lowres_cortex.faces lowres_cortex.normals] = load_tobj(outputFile);
                            lowres_cortex.facevertexcdata = zeros(size(lowres_cortex.vertices,1),1);

                            fprintf('  Saving\n');
                            eval(['save ' fileName ' lowres_cortex']);
                            clear lowres_cortex;
                    end
            end
        end
    end
end

%%  Find values for projection
%   This function is called in ProjectElectrodes() to find the minimum
%   distance on a surface "gs" that is closest for each electrode in "els"

function out_els = p_zoom(els, gs, index, checkdistance, ignoreTrodes)
    %function out_els=p_zoom(els, gs, index);
    % this function finds the minimum distance on a surface "gs" (pts x 3) that is
    % closest for each electrode in "els" (N x 3)
    % index indicates the number of point for calculation of norm, 0 if global
    % checkdistance = 1 to indicate that electrodes within 3 mm of the surface
    % are not projected

    %use a local set of electrodes to determine the orthogonal direction of a
    %given electrode? -- enter "0" for global, and number to use otherwise
    %   Created by:
    %   D. Hermes & K.J. Miller 
    %   Dept of Neurology and Neurosurgery, University Medical Center Utrecht
    %
    %   Version 1.1.0, released 26-11-2009

    if ~exist('ignoreTrodes','var')
        ignoreTrodes = [];
    end

    lcl_num = index;
    checkdistance_dist = 3; % 3 mm

    if checkdistance == 2
        disp('electrodes projected to closest point, no norm')
    end

    if mean(els(:,1))<0
        disp('left grid');
        % delete right hemisphere
        % gs=gs(gs(:,1)<=0,:);
    else
        disp('right grid');
        % delete right hemisphere
        % gs=gs(gs(:,1)>=0,:);
    end

    if lcl_num == 0 %global estimate of principal direction most orthogonal to array
        [v,d] = eig(cov(els)); %all vecs
        nm = v(:,find(diag(d) == min(diag(d)))); %vec we want
        nm = nm*sign(nm(1)*mean(els(:,1)));%check for left or rigth brain, invert nm if needed
    end

    out_ind = zeros(size(els,1),1);
    out_ind_rev = zeros(size(els,1),1);
    for k = 1:size(els,1)

        if find(ignoreTrodes == k,1,'first')>0
            continue;
        end
        
        %sub array?
        if lcl_num > 0, % get principal direction most orthogonal to sub-array
            [y,ind] = sort(dist(els,els(k,:)'),'ascend');%select closest for sub-array
            [v,d] = eig(cov(els(ind(1:lcl_num),:))); %all vecs
            nm = v(:,find(diag(d)==min(diag(d)))); %vec we want
            nm = nm*sign(nm(1)*mean(els(:,1)));%check for left or right brain, invert nm if needed
        end
        
        npls = [gs(:,1)-els(k,1) gs(:,2)-els(k,2) gs(:,3)-els(k,3)]; %x,y,z lengths
        % calculate distance
        npls_dist = sqrt(sum((npls).^2,2));
        % check whether distance is < 3 mm
        distancesm2 = 0;
        if npls_dist(find(npls_dist == min(npls_dist)),:) < checkdistance_dist
            %disp(['distance < 3 mm electrode ' int2str(k) ]);
            distancesm2 = 1;
        end

        if checkdistance == 1 && distancesm2 == 1 % electrode too close to surface to project
            out_ind(k) = find(npls_dist == min(npls_dist),1); %find minimum distance
        elseif checkdistance == 2
            out_ind(k) = find(npls_dist == min(npls_dist),1); %find minimum distance
        else
            npls_unit = npls./repmat((sum(npls.^2,2).^.5),1,3); % normalize npls to get unit vector
            npdot = (npls_unit*nm); %distance along eigenvector direction (dot product)
            % only take gs within 2.5 cm distance
            npdot(npls_dist>25) = 0;
            %npdotrev=(npls_unit*-nm); % distance in reverse eigenvector direction
            [a b] = find(abs(npdot) == max(max(abs(npdot))),1);
            out_ind(k) = a;%find minimum distance, max dot product
            %out_ind_rev(k)=find(npdotrev==max(npdotrev),1); %find minimum distance, max dot product
        end
    end
    out_ind(out_ind == 0) = [];
    out_els = els;
    out_els(setdiff(1:size(els,1),ignoreTrodes),:) = gs(out_ind,:);
    % out_els_rev=gs(out_ind_rev,:);
    % % check distance to new electrodes
    % out_els_dist=sqrt(sum((els-out_els).^2,2));
    % out_els_rev_dist=sqrt(sum((els-out_els_rev).^2,2));

    % plot on surface to check
    figure
    plot3(els(:,1),els(:,2),els(:,3),'r.','MarkerSize',20);
    hold on;
    plot3(out_els(:,1),out_els(:,2),out_els(:,3),'g.','MarkerSize',20);
    plot3(gs(:,1),gs(:,2),gs(:,3),'k.','MarkerSize',1);
    axis equal;

    for i=1:size(els,1);
        el = els(i,:);
        text(el(1),el(2),el(3),sprintf('%i',i));
    end
end

%%  Project electrodes
%   This function projects the electrodes onto the generated surface

function ProjectElectrodes(subjDir )
    BioImageToMatlab(getenv('recon_patientCode'));

    load(fullfile(subjDir, 'bis_trodes.mat'));

    hemi = 'r';

    originOffset = [0 0 0];

    for arrayName = TrodeNames
        eval(sprintf('target = %s;', arrayName{:}));
        
        hemi = [];
        fprintf('Processing ''%s''...\n',arrayName{:});
        while isempty(hemi)
            fprintf('  Hemisphere (r/l): \n');
            choice = input('    => ','s');
            switch choice
            case 'r'
                hemi = 'rh';
            case 'l'
                hemi = 'lh';
            otherwise
                error('  Bad choice\n');
            end
        end

        fprintf('  Type of grid:\n',arrayName{:});
        fprintf('    %5i - M X N - Grid\n',1);
%         fprintf('%5i - 2 x N - Wide Strip\n',2);

        fprintf('    %5i - 1 X N - Thin Strip\n',2);
        fprintf('    %5i - 1 X N - Depth\n',3);
        choice = input('  => ');

        switch choice
            case 1
                index = 5;
                checkDistance = 1;
            case 2
                index = 0;
                checkDistance = 2;
            case 3
                continue;
            otherwise
                error('  Bad choice\n');
        end

        load([subjDir '/surf/' getenv('recon_patientCode') '_cortex_' hemi '_hires.mat']);

        [x,y,z]=ind2sub(size(brainHull),find(brainHull>0)); 
        % from indices 2 native
        gs=([x y z]*transformMat(1:3,1:3)')+repmat(transformMat(1:3,4),1,length(x))';
        
        out_els = [];
        ignoreEls = [];

        while(1)
            out_els=p_zoom(target,gs,index,checkDistance, ignoreEls);
    %         out_els(ignoreEls,:) = target(ignoreEls,:);
    %         out_els=p_zoom(target + repmat(originOffset,size(target,1),1),gs,index,checkDistance);

            if hemi=='l'
                view(240, 30);     
            elseif hemi=='r'
                view(60, 30);      
            end

            % errors = sqrt((out_els(:,1) - Grid(:,1)).^2 + (out_els(:,2) - Grid(:,2)).^2 + (out_els(:,3) - Grid(:,3)).^2);
            % bad = find((errors - mean(errors)) ./ std(errors) > 2);

        %     ctmr_gauss_plot(cortex,[0 0 0],0,hemi); 
        %     hold on;
        %     plot3(target(:,1),target(:,2),target(:,3),'markersize',20,'linestyle','none','marker','.','color','r');
        %     plot3(out_els(:,1),out_els(:,2),out_els(:,3),'markersize',20,'linestyle','none','marker','.','color','g');

            fprintf('Options:\n');
            fprintf('  1. Projected positions\n');
            fprintf('  2. Original Positions\n');
            fprintf('  3. Smooth electrodes and re-project\n');
            fprintf('  4. Change origin offset\n');
            fprintf('  5. Ignore projection for subset of eletrodes\n');
            fprintf('  6. Manually set electrode position\n');
            choice = lower(input('Default [1] => ','s'));

            if isempty(choice)
                choice = '1';
            end

            switch choice
                case '1'
                    eval(sprintf('%s = out_els;', arrayName{:}));
                    break;
                case '2'
                    fprintf('Ignoring projection.\n');
                    break;
                case '4'
                    originOffset = input('Change origin offset: \','s');
                    originOffset = str2num(originOffset);
                case '3'
                    fprintf('Smoothing electrodes...\n');

                    temp = target;
                    target = zeros(size(temp));
                    
                    if (checkDistance == 1) % working with a grid, % THIS ASSUMES 8x8!!
                        target(1,:) = temp(1,:); target(8,:) = temp(8,:); target(57,:) = temp(57,:); target(64,:) = temp(64,:);

                        % smooth edges
                        for i=2:7; target(i,:) = (temp(i-1,:) + temp(i+1,:)) ./ 2; end
                        for i=58:63; target(i,:) = (temp(i-1,:) + temp(i+1,:)) ./ 2; end
                        for i=2:7; target(i*8,:) = (temp((i-1)*8,:) + temp((i+1)*8,:)) ./ 2; end
                        for i=2:7; target(i*8-7,:) = (temp((i-1)*8-7,:) + temp((i+1)*8-7,:)) ./ 2; end

                        % smooth middle
                        for row=2:7
                            for col=2:7
                                target(((row-1)*8)+col,:) = (temp(((row-2)*8)+col,:) + temp(((row)*8)+col,:) + temp(((row-1)*8)+col+1,:) + temp(((row-1)*8)+col-1,:)) ./ 4;
                            end
                        end
                        
                    elseif (checkDistance == 2) % working with a strip
                        target(1,:) = temp(1,:);
                        target(end,:) = temp(end,:);
                        
                        for c = 2:(size(target,1)-1)
                            target(c,:) = (temp(c-1,:) + temp(c+1,:)) ./ 2; 
                        end
                    else
                        fprintf(' .. couldn''t smooth\n');
                    end
                    
                case '5'
                    fprintf('Add electrodes to be ignored in vector format.  i.e. [1 2 4] or [1:5 8:11] etc\n');
                    elecsChosen = input('=> ','s');
                    eval(sprintf('ignoreEls = unique(union(%s, ignoreEls));',elecsChosen));
                case '6'
                    fprintf('Enter electrode index: \n');
                    trodeIdx = input('=>', 's');
                    trodeIdx = str2num(trodeIdx);
                    
                    fprintf('Enter new coordinates i.e. [20 -14.3 7.8]: \n');
                    coords = input('=>', 's');
                    eval(sprintf('target(%d,:) = [%s];', trodeIdx, coords));
                otherwise
                    fprintf('Bad choice\n');
            end
            close all;
        end
        close all;
    end

    AllTrodes = [];
    for name = TrodeNames
        name = name{1};
        eval(sprintf('AllTrodes = cat(1,AllTrodes,%s);', name));
    end

    cmd = ['save ' subjDir 'trodes.mat AllTrodes TrodeNames '];
    for name = TrodeNames
        name = name{1};
        cmd = [cmd name ' '];
    end

    cmd = [cmd ';'];
    eval(cmd);
end

