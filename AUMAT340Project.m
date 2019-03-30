%{
Philippe Nadon

Project
AUMAT 340
March 12, 2019
%}

function [] = AUMAT340Project()

    % Choose parts to run corresponding Matlab code
    parts = [1, 2, 3, 4, 5, 6];
    year = 2100;
    
    expGuess = [-0.4, 0.02, 0.1];
    expFun = @(x, t) x(1) + x(2) * exp( x(3) * t);
    
    mixGuess1 = [-0.5, 0.02, -0.002, -0.4];
    mixFun1 = @(x, t) x(1) + (x(3) * t + x(4)) .* exp( x(2) * t);
    
    mixGuess2 = [-0.5, 0.02, 0.03, -0.002, -0.4];
    mixFun2 = @(x, t) x(1) + (x(3) * t .^ 2 + x(4) * t + x(5)) .* exp( x(2) * t);

    file = readtable("data.txt", 'HeaderLines', 77);
    data = table2array( file(:,5) );
    N = length(data);
    disp(N);
    t = zeros( N, 1);
    for i = 1:N
        t(i) = floor( (i + 4) / 12);
    end
    options = optimoptions('lsqcurvefit', 'MaxFunctionEvaluations', 40000, 'MaxIterations', 40000);
    
    for part = parts
        disp(".");
        switch part
            case 1
                disp("|======= Linear Result =====================|");
                res = extrapolatePolyToYear( t, data, 1);
                displayPolyRes( t, data, res, year);
            case 2
                disp("|======= Quadratic Result ==================|");
                res = extrapolatePolyToYear( t, data, 2);
                displayPolyRes( t, data, res, year);
            case 3
                disp("|======= Cubic Result ======================|");
                res = extrapolatePolyToYear( t, data, 3);
                displayPolyRes( t, data, res, year);
            case 4
                disp("|======= Exponential Result ================|");
                res = lsqcurvefit( expFun, expGuess, t, data, [], [], options);
                displayRes( expFun, t, data, res, year);
            case 5
                disp("|=== Exponential and Degree-1 Polynomial ===|");
                res = lsqcurvefit( mixFun1, mixGuess1, t, data, [], [], options);
                displayRes( mixFun1, t, data, res, year);
            case 6
                disp("|=== Exponential and Degree-2 Polynomial ===|");
                res = lsqcurvefit( mixFun2, mixGuess2, t, data, [], [], options);
                displayRes( mixFun2, t, data, res, year);
        end
    end
end

function res = extrapolatePolyToYear( t, data, degree)
    N = length(data);
    rows = degree + 1;
    A = zeros( rows, rows);
    B = zeros( rows, 1);
    
    for n = 1:N
        for i = 1:rows
            % for A:
            % The diagonal entries:
            A(i, i) = A(i, i) + t(n) ^ (2 * ( i - 1));
            % The upper / lower entries (mirror along diagonal)
            for j = 1:(i - 1)
                A( i, j) = A(i, j) + t(n) ^ ( i + j - 2);
                A( j, i) = A(i, j);
            end
            
            % for B:
            B(i) = B(i) + data(n) * t(n) ^ (i - 1);
        end
    end
    disp(A);
    disp(B);

    res = linsolve( A, B);
end

function [] = displayPolyRes( t, data, res, year)
    switch( length(res))
        case 1
            fprintf("alpha: %d\n", res(1));
        case 2
            fprintf("alpha: %d, beta: %d\n", res(1), res(2));
        case 3
            fprintf("alpha: %d, beta: %d, gamma: %d\n", res(1), res(2), res(3));
        case 4
             fprintf("alpha: %d, beta: %d, gamma: %d, delta: %d\n", res(1), res(2), res(3), res(4));
    end
    
    numPoints = length(t) + year - 1850 - t( length(t));
    tExtended = zeros( numPoints, 1);
    resPlot = zeros( numPoints, 1);
    
    for i = 1:numPoints
        if i <= length(t)
            tExtended(i) = t(i);
        else 
            tExtended(i) = tExtended(i - 1) + 1;
        end
        for j = 1:length(res)
            resPlot(i) = resPlot(i) + res(j) * tExtended(i) ^ (j - 1); 
        end
    end
    
    errorSum = 0;
    for i = 1:length(data)
    errorSum = errorSum + (data(i) - resPlot(i))^2;
    end
    
    fprintf("Anomaly at year %d: %d\n", year, resPlot( numPoints));
    fprintf("Cumulative Error: %d\n", errorSum);
    
    t = arrayfun( @(x) x + 1850, t);
    tExtended = arrayfun( @(x) x + 1850, tExtended);
    
    figure
    plot( t, data);
    %set(gca,'FontSize',24)
    xlabel('Year');
    ylabel('Anomaly');
    hold on
    plot( tExtended, resPlot);
    hold off
    
end

function [] = displayRes( fun, t, data, res, year)
    disp(res);
    
    numPoints = length(t) + year - 1850 - t( length(t));
    tExtended = zeros( numPoints, 1);
    resPlot = zeros( numPoints, 1);
    
    for i = 1:numPoints
        if i <= length(t)
            tExtended(i) = t(i);
        else 
            tExtended(i) = tExtended(i - 1) + 1;
        end
        resPlot(i) = fun( res, tExtended(i));
    end
    
    errorSum = 0;
    for i = 1:length(data)
    errorSum = errorSum + (data(i) - resPlot(i))^2;
    end
    
    fprintf("Anomaly at year %d: %d\n", year, resPlot( numPoints));
    fprintf("Cumulative Error: %d", errorSum);
    
    t = arrayfun( @(x) x + 1850, t);
    tExtended = arrayfun( @(x) x + 1850, tExtended);
    
    figure
    plot( t, data);
    %set(gca,'FontSize',24)
    xlabel('Year');
    ylabel('Anomaly');
    hold on
    plot( tExtended, resPlot);
    hold off
end