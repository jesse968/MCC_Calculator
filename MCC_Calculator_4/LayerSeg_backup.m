function [ImSeg, Layer] = LayerSeg( Image, Layer )
    
    Layer = getRetinalLayers( Image, Layer );
    
    % Convert grayscale image to RGB
    rgbImage = cat(3, Image, Image, Image);

    ImSeg = insertShape( rgbImage, 'line', [Layer(:,1) Layer(:,2)], ...
               'Color', [255, 0, 0], 'LineWidth', 1, 'Opacity', 1);
end

function Layer = getRetinalLayers( img, Layer )
        
    % constants used for defining the region for segmentation of individual layer
    dark2bright   = true;
    regionOffset  = 10;
    topOffset     = 10;
    topRightExtension = 0.4;
    
    % smooth image with specified kernels, for denosing
    img = imfilter(img, fspecial('gaussian', [5 10], 5), 'replicate');
    
    roiImg = GetRoiImg(img, Layer, regionOffset, topOffset, topRightExtension);

    % create adjacency matrices and its elements base on the image.
    matrixWeight = getAdjacencyMatrix(img, roiImg, dark2bright);
    
    [layer, ~] = getRetinalLayersCore( img, matrixWeight, roiImg );

    Layer = layer;  % 将当前图像的分层结果用于下一组图像的基准值
end

function MatrixWeight = getAdjacencyMatrix(inputImg, RoiImg, Dark2Bright)
    
    szImg = size(inputImg);
    
    % get vertical gradient image
    [~, gradImg] = gradient(double(inputImg),1,1);
    gradImg = -1*gradImg;
    
    % normalize gradient
    gradImg = (gradImg-min(gradImg(:)))/(max(gradImg(:))-min(gradImg(:)));
    
    minWeight = 1E-5;

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
    if true == Dark2Bright         % Default
        gradImg = gradImg*-1+1;    % get the "invert" of the gradient image
    end

    adjWeight = 2 - gradImg(roiIndices(:)) - gradImg(roiIndicesAdj(:)) + minWeight;
    
    % pad minWeight on the side
    imgTmp = nan(size(gradImg));
    imgTmp(:,1)   = 1;
    imgTmp(:,end) = 1;
    idx = int32(size(gradImg,2) * 0.6);
    imgTmp(11, idx:end) = 1;

    imageSideInd = ismember(roiIndicesAdj,find(imgTmp(:)==1));
    adjWeight(imageSideInd) = minWeight;

    MatrixWeight(:,1) = roiIndices;
    MatrixWeight(:,2) = roiIndicesAdj;
    MatrixWeight(:,3) = adjWeight;
end

function [rPaths, Img] = getRetinalLayersCore(Img, MatrixWeight, RoiImg)

    adjMA  = MatrixWeight(:,1);
    adjMB  = MatrixWeight(:,2);
    adjMmW = MatrixWeight(:,3);
    
    % include only region of interst in the adjacency matrix
    includeA = ismember(adjMA, find(RoiImg(:) == 1));
    includeB = ismember(adjMB, find(RoiImg(:) == 1));
    keepInd = includeA & includeB;
    keepInd = ~keepInd;
    
    %get the shortestpath
    adjMmW(keepInd) = 30;
    adjMatrixMW = sparse( adjMA, adjMB, adjMmW, numel(Img(:)), numel(Img(:)) );
    G = digraph(adjMatrixMW);

    [ path, ~ ] = shortestpath( G, 1, numel(Img(:)) );
    
    % convert path indices to subscript
    [rPathsNew(:,2), rPathsNew(:,1)] = ind2sub( size(Img), path );
    
    % save data
    rPaths = rPathsNew;
end

function RoiImg = GetRoiImg(Img, Layer, regionOffset, TopOffset, TopRightExtension)
    % initialize region of interest
    szImg  = size(Img);
    RoiImg = zeros(szImg);
    
    % avoid the top part of image
    RoiImg(1:TopOffset,:) = 0;

    % 快速通道
    RoiImg(:, 1)   = 1;
    RoiImg(:, end) = 1;
    index = uint32(szImg(2) * (1-TopRightExtension));
    RoiImg(TopOffset+1, index:end) = 1;
    
    x = Layer(:,1);
    y = Layer(:,2);

    % select region of interest based on reference layer
    prvIdx = mean(y(x == 2));
    for k = 1:szImg(2)
        index = mean(y(x == k));
        diff  = abs(index-prvIdx);
        if( diff > (2*regionOffset+1))
            regionOffset = ceil(diff / 2);
        end
        prvIdx = index;

        startInd = uint32(index - regionOffset);
        endInd   = uint32(index + regionOffset);

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
        RoiImg(startInd:endInd,k) = 1;
    end
end