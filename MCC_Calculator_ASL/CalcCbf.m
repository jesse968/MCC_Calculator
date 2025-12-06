function [Result, Graph] = CalcCbf( Slices, Reduction, FrameRate )
    
    frames = numel(Slices);
    [height, width, ~] = size(Slices{1});

    images = zeros( frames, height, width );
    for i=1:frames
        images(i,:,:) = rgb2gray(Slices{i});
    end

    cbf = zeros(frames,width);
    for j=1:height
        K(:,:) = images(:,j,:);
        frq = fft(double(K(:,:)), frames);
        cbf=cbf+abs(frq);
    end

    s = sum(cbf');

    index = round(Reduction*frames);
    s(1:index) = min(s(index:end));
    
    [mx1, mx2] = max(s(1:round(frames/2)));

    coeff = round(frames / FrameRate);
    Result = mx2/coeff;
    Graph  = s(1:round(frames/2));
end
