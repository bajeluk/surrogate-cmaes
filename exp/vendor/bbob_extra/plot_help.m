function [] = plot_help()
    handles = benchmarksextra('handles');
    testFn = handles(5);
    testFn = cell2mat(testFn);
    testFn('init', [], 355);
    func = @(x, y) testFn([x; y; +1.9e+00; -4.8e+00]);
    
    res = [];
    x = -5:0.1:5;
    y = -5:0.1:5;
    for i = y
        row = [];
        for j = x
            row = [row, func(i, j)];
        end
        res = [res; row];
    end
    display(res);
    
    f = figure('visible', 'off');
    mesh(x, y, res)
    print -djpeg ../prubehy_fci/f205-samp05-355/3-4-fixed-mesh.jpg
    
    contour(x, y, res)
    print -djpeg ../prubehy_fci/f205-samp05-355/3-4-fixed-contour.jpg
    close(f)
   
    % contour(x, y, res)
    
end