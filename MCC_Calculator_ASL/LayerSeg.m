function [ImSeg, Layer] = LayerSeg( Image, Layer, AuxLineDist, SegOffset )
    
    Layer = getRetinalLayers( Image, Layer, SegOffset );
    [height,~] = size(Image);
    
    % Convert grayscale image to RGB
    rgbImage = cat(3, Image, Image, Image);

    smoothSignal = true;
    if true == smoothSignal
        % 尝试平均滤波
        windowSize = 5;
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        smoothedY = filter(b, a, Layer(:,2)); % 平滑Layer的纵坐标
        Layer(:,2) = smoothedY;
    end

    ImSeg = insertShape( rgbImage, 'line', [Layer(:,1) Layer(:,2)], ...
               'Color', [255, 0, 0], 'LineWidth', 2, 'Opacity', 1);

    if 0 ~= AuxLineDist
        newLayer(:,1) = Layer(:,1);
        newLayer(:,2) = min( height, Layer(:,2)+AuxLineDist );
        newLayer(:,2) = max( 1, Layer(:,2)+AuxLineDist );
    
        ImSeg = insertShape( ImSeg, 'line', [newLayer(:,1) newLayer(:,2)], ...
                   'Color', [0, 255, 0], 'LineWidth', 1, 'Opacity', 1);
    end
end

function Layer = getRetinalLayers( Img, Layer, SegOffset )
        
    % constants used for defining the region for segmentation of individual layer
    dark2bright   = true;
    topOffset     = 10;

    topExtension = [0.08 0];
    
    % smooth image with specified kernels, for denosing
    % img = imfilter(img, fspecial('gaussian', [5 10], 5), 'replicate');
    filterSize = [1, 3];
    Img = medfilt2(Img, filterSize);
    
    roiImg = GetRoiImg(Img, Layer, SegOffset, topOffset, topExtension);

    % create adjacency matrices and its elements base on the image.
    matrixWeight = getAdjacencyMatrix(Img, roiImg, dark2bright, topExtension);
    
    [layer, ~] = getLayersCore( Img, matrixWeight );

    if ~isempty(layer)
        Layer = layer;  % 将当前图像的分层结果用于下一组图像的基准值
    end
end

function MatrixWeight = getAdjacencyMatrix(inputImg, RoiImg, Dark2Bright, TopExtension)
    
    szImg = size(inputImg);
    
    % get vertical gradient image

    inputImg = double(inputImg);
    gradImg  = inputImg(3:end,:) - inputImg(1:end-2,:);
    gradImg  = vertcat(zeros(1,szImg(2)), gradImg, zeros(1,szImg(2)));

    % generate adjacency matrix, see equation 1 in the refered article.
    neighborIterX = [1 1  1 0  0 -1 -1 -1];
    neighborIterY = [1 0 -1 1 -1  1  0 -1];
    
    % 提取roi image中的坐标索引
    roiIndices = find(RoiImg==1);

    % convert adjMA to subscripts
    [adjMAx,adjMAy] = ind2sub(szImg,roiIndices);

    roiIndices = roiIndices';
    roiIndSize = size(roiIndices);
    
    % prepare to obtain the 8-connected neighbors of adjMAsub
    % repmat to [1,8]
    neighborIterX = repmat(neighborIterX, [roiIndSize(1),1]);
    neighborIterY = repmat(neighborIterY, [roiIndSize(1),1]);
    
    % repmat to [8,1]
    roiIndices = repmat(roiIndices,[1 8]);
    adjMAx = repmat(adjMAx, [1 8]);
    adjMAy = repmat(adjMAy, [1 8]);
    
    % get 8-connected neighbors of adjMAsub
    % adjMBx,adjMBy and adjMBsub
    adjMBx = adjMAx+neighborIterX(:)';
    adjMBy = adjMAy+neighborIterY(:)';
    
    % make sure all locations are within the image.
    keepInd = adjMBx > 0 & adjMBx <= szImg(1) & ...
              adjMBy > 0 & adjMBy <= szImg(2);

    roiIndices = roiIndices(keepInd);
    adjMBx = adjMBx(keepInd);
    adjMBy = adjMBy(keepInd);
    
    roiIndicesAdj = sub2ind(szImg,adjMBx(:),adjMBy(:))';
    
    % calculate weight
    minWeight = 0.001;
    if true == Dark2Bright
        gradImg(gradImg<0) = 0;
        adjWeight = gradImg(roiIndices(:)) + gradImg(roiIndicesAdj(:));
        adjWeight = normalize(adjWeight, 'range', [0, 1]);
        gamma = 2;
        adjWeight = adjWeight.^gamma;
        adjWeight = abs(adjWeight - 1);
        
        % 放大比重在特定范围内
        % 将介于0.9到1之间的数值放大20倍
        condition = (adjWeight >= 0.9) & (adjWeight <= 1);
        adjWeight(condition) = adjWeight(condition) * 100;

        % 将介于0.8到0.9之间的数值放大5倍
        condition = (adjWeight >= 0.8) & (adjWeight <= 0.9);
        adjWeight(condition) = adjWeight(condition) * 5;

        adjWeight = adjWeight + minWeight;
        %gradImg = normalize(gradImg, 'range', [0, 1]);
        %gamma = 2.0;
        %gradImg = gradImg.^gamma;
        %gradImg = abs(gradImg - 1);
    end

    % adjWeight = gradImg(roiIndices(:)) + gradImg(roiIndicesAdj(:));
    
    % pad minWeight on the side (暂时关掉两边的快速通道)
    imgTmp = nan(size(gradImg));
    imgTmp(:,1)   = 1;
    imgTmp(:,end) = 1;

    leftExt  = int32(size(gradImg, 2) * TopExtension(1));
    rightExt = int32(size(gradImg,2 ) * (1-TopExtension(2)));
    
    imgTmp(11, 1:leftExt) = 1;
    imgTmp(11, rightExt:end) = 1;

    imageSideInd = ismember(roiIndicesAdj,find(imgTmp(:)==1));
    adjWeight(imageSideInd) = minWeight;

    MatrixWeight(:,1) = roiIndices;
    MatrixWeight(:,2) = roiIndicesAdj;
    MatrixWeight(:,3) = adjWeight;
end

function [rPaths, Img] = getLayersCore(Img, MatrixWeight)

    adjMA  = MatrixWeight(:,1);
    adjMB  = MatrixWeight(:,2);
    adjMmW = MatrixWeight(:,3);
    
    % include only region of interst in the adjacency matrix
    %includeA = ismember(adjMA, find(RoiImg(:) == 1));
    %includeB = ismember(adjMB, find(RoiImg(:) == 1));
    %keepInd = includeA & includeB;
    %keepInd = ~keepInd;
    
    %get the shortestpath
    %adjMmW(keepInd) = 30;
    adjMatrixMW = sparse(adjMA, adjMB, adjMmW, numel(Img(:)), numel(Img(:)));
    G = digraph(adjMatrixMW);

    [path, ~] = shortestpath(G, 1, numel(Img(:)), 'Method', 'positive');
    
    if isempty(path)
        sprintf("Error");
        return;
    end

    % convert path indices to subscript
    [rPathsNew(:,2), rPathsNew(:,1)] = ind2sub(size(Img), path);
    
    % save data
    rPaths = rPathsNew;
end

function RoiImg = GetRoiImg(Img, Layer, regionOffset, TopOffset, TopExtension)
    % initialize region of interest
    szImg  = size(Img);
    RoiImg = zeros(szImg);
    
    % avoid the top part of image
    RoiImg(1:TopOffset,:) = 0;

    % 快速通道
    RoiImg(:, 1)   = 1;
    RoiImg(:, end) = 1;

    leftExt  = uint32(szImg(2) * TopExtension(1));
    rightExt = uint32(szImg(2) * (1-TopExtension(2)));

    RoiImg(TopOffset+1:TopOffset+6, 1:leftExt)    = 1;
    RoiImg(TopOffset+1:TopOffset+6, rightExt:end) = 1;
    
    x = Layer(:,1);
    y = Layer(:,2);

    % select region of interest based on reference layer
    prvIdx = mean(y(x == 2));
    for k = 1:size(x,1)
        index = y(k);
        diff  = abs(index-prvIdx);
        if( diff > (2*regionOffset+1))
            offset = ceil(diff / 2);
        else
            offset = regionOffset;
        end
        prvIdx = index;

        startInd = uint32(index - offset);
        endInd   = uint32(index + offset);

        % error checking
        if startInd > endInd
            startInd = endInd - 1;
        end            
        
        if startInd <= TopOffset
            startInd = TopOffset + 1;
        end
        
        if endInd > szImg(1)
            endInd = szImg(1);
        end
                        
        % set region of interest at column k from startInd to endInd
        leftIdx  = max(1, x(k)-3);
        rightIdx = min(szImg(2), x(k)+3);
        RoiImg(startInd:endInd,leftIdx:rightIdx) = 1;
    end
end