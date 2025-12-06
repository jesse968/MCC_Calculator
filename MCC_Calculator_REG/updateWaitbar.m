function updateWaitbar(h, totalIterations, reset)
    persistent count

    if true == reset
        count = 0;
        return;
    end

    if isempty(count)
        count = 0;
    end
    count = count + 1;
    progress = count / totalIterations;
    waitbar(progress, h);
end