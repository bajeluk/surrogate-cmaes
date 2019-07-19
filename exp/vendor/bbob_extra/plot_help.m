function [] = plot_help()
    handles = benchmarksextra('handles');
    testFn = handles(6);
    testFn = cell2mat(testFn);
    testFn('init', [], 1);
    func = @(x, y) testFn([x; y]);
    
    res = [];
    x = -5:0.02:5;
    y = -5:0.04:5;
    for i = y
        row = [];
        for j = x
            row = [row, func(j, i)];
        end
        res = [res; row];
    end
    display(res);
    
    f = figure('visible', 'off');
    mesh(x, y, res)
    print -djpeg ~/School/diplomka-support/prubehy_fci/f206-detailed/resize-1-mesh.jpg
    
    contour(x, y, res)
    print -djpeg ~/School/diplomka-support/prubehy_fci/f206-detailed/resize-1-contour.jpg
    close(f)
   
    % contour(x, y, res)
    
end