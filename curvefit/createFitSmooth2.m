function [fitresult, gof] = createFitSmooth2(x, y, w, smooth, outliers)
%  Create a fit.
%      X Input : x
%      Y Output: y
%      Weights : w
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%--------------------------------------------------------------------------
[xData, yData, weights] = prepareCurveData( x, y, w );

% Set up fittype and options.
%fo = fitoptions('SmoothingParam',smooth, 'Weights',weights, 'Exclude',outliers);
%ft = fittype('smoothingspline', 'options',fo);

%opts. = smooth; %0.00632767277496633;
%opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, 'smoothingspline', 'SmoothingParam',smooth, 'Weights',weights, 'Exclude',outliers);

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'y vs. x with w', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'x' );
% ylabel( 'y' );
% grid on


