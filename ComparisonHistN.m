function [hFig] = ComparisonHistN(cvfData, vfBins, cfPips, vfP, cvfColors, strMeasureName, vbFillBars)

% ComparisonHistN - FUNCTION Make a histogram comparing an arbitrary number of data sets
%
% Usage: [hFig] = ComparisonHist3({vfA, vfB, vfC, ...}, vfBins, ...
%                                 <{fPipA, fPipB, fPipC, ...}, ...>
%                                 <vfP> or <fhTest> or <{fhTest, nOutputInd}>
%                                 <{vfColorA, vfColorB, vfColorC, ...}, ...>
%                                 <strMeasureName, vbFillBars>)
%
% 'vfA', 'vfB', ... are vectors of data, which will be plotted as
% proportional histograms on a single set of axes. 'vfBins' is the list of
% bin *EDGES*, as used in histc.
%
% The optional arguments 'fPipA', 'fPipB', ... are scalar values
% corresponding to each vector of data, and will be plotted as "pips" above
% the histograms. If they are not supplied, pips will be plotted as medians
% of each vector. If an empty matrix is supplied, that pip will not be
% displayed.
%
% The optional argument 'vfP' is a vector of P values indicating the result
% of a statistical test of some pairwise difference between two of the data
% vectors. It has the format ['fP_AB' 'fP_AC' 'fP_BC' ...]. Significant
% differences will be indicated by asterisks and bars joining the pips of
% the appropriate distributions.  If 'vfP' is not provided, then the
% Wilcoxon rank-sum test of differences in medians will be used.
%
% Alternatively, a function handle can be provided in 'fhTest'. This
% function must return the P value indicating the significance of a
% difference between two distributions: fP = @(vfX, vfY). This function
% will be applied pair-wise between each data vector. NOTE: ttest2 does NOT
% return the P value as the first argument, and so cannot be used directly
% in this way. In this case, pass the function and which output to take as
% a cell array: {@ttest2, 2}. This form will take the second output of
% ttest2 as the P value to use.
%
% The optional arguments 'vfColorA', 'vfColorB', ... set the
% colours for the corresponding histogram. If not supplied, a set of
% colours will be chosen for you. The first data set will be plotted in
% black.
%
% The optional argument 'strMeasureName' is used as the x label of the
% figure.
%
% The optional argument 'vbFillBars' is a boolean vector, with each element
% corresponding to a data set 'vfA', 'vfB', ...  A value of 'true'
% indicates that the corresponding data set should be plotted using filled
% histogram bars. A value of 'false' indicates that the corresponding data
% set should be plotted using a stair plot. By default, a stair plot is
% used for all data sets except the first.
%
% 'hFig' will be the figure handle of the comparison figure. The title of
% the figure will contain useful statistical information about the
% histograms.

% Author: Dylan Muir <dylan.muir@unibas.ch>
% Created: 2013

% -- Default arguments

DEF_strMeasureName = 'Measure';
DEF_fP = @(a,b)ranksum(a(~isnan(a)), b(~isnan(b)));


%% -- Clean and reshape input data

if (~iscell(cvfData))
   cvfData = {cvfData};
end
cvfData = cvfData(~cellfun(@isempty, cvfData));
cvfData = cellfun(@(c)reshape(c, [], 1), cvfData, 'UniformOutput', false);


%% -- Check arguments

if (~exist('cfPips', 'var'))
   cfPips = cellfun(@nanmedian, cvfData, 'UniformOutput', false);
end

if (numel(cfPips) ~= numel(cvfData))
   error('*** ComparisonHistN: Error: Number of pips does not match number of data sets.');
end

if (~exist('vfP', 'var'))
   vfP = DEF_fP;
end

% - Convert cell vfP argument to function handle
if iscell(vfP)
   vfP = @(varargin)(outputs(vfP{1}, vfP{2}, varargin{:}));
end

% - Evaluate P values using supplied function
if (isa(vfP, 'function_handle'))
   nPEntry = 1;
   vfThesePs = [];
   for (nSeries1 = 1:numel(cvfData))
      for (nSeries2 = nSeries1+1:numel(cvfData))
         vfThesePs(nPEntry) = vfP(cvfData{nSeries1}, cvfData{nSeries2}); %#ok<AGROW>
         nPEntry = nPEntry + 1;
      end
   end
   vfP = vfThesePs;
end

% - Get a set of colours to use for plots
if (~exist('cvfColors', 'var') || isempty(cvfColors))
   nNumColours = numel(cvfData);
   
   if (nNumColours > 1)
      if (exist('pmkmp', 'file'))
         cvfColors = pmkmp(nNumColours, 'IsoL');
      else
         cvfColors = hsv(nNumColours);
      end
         
      cvfColors = mat2cell(cvfColors, ones(1, nNumColours));
   end
   
   % - First series is black
   cvfColors{1} = [0 0 0];
end

% - Set a default measure name
if (~exist('strMeasureName', 'var'))
   strMeasureName = DEF_strMeasureName;
end

% - By default, don't fill histogram bars for 2nd and subsequent series
if (~exist('vbFillBars', 'var') || isempty(vbFillBars))
   vbFillBars = false(size(cvfData));
   vbFillBars(1) = true;
end


%% -- Compute histograms

% - Compute histogram counts for each data set
for (nSeries = numel(cvfData):-1:1)
   cvfHist{nSeries} = histc(cvfData{nSeries}, vfBins) ./ nnz(~isnan(cvfData{nSeries})) * 100;
   cvfHist{nSeries}(end-1) = sum(cvfHist{nSeries}(end-1:end));
   cvfHist{nSeries} = cvfHist{nSeries}(1:end-1);
end

% - Work out level at which pips should be plotted
fPipLevel = nanmax(cellfun(@nanmax, cvfHist)) * 1.1;

% - Check significance
cstrSig = arrayfun(@check_sig, vfP, 'UniformOutput', false);


%% -- Plot histograms

% - Make a figure, if necessary
newplot;
bHold = ishold;
hFig = gcf;

% - Plot histogram data sets
for (nHist = 1:numel(cvfData))
   if (vbFillBars(nHist))
      % - Plot a filled series
      vhBar(nHist) = bar(vfBins, [cvfHist{nHist}(:)' nan], 'histc'); %#ok<AGROW>
      set(vhBar(nHist), 'EdgeColor', cvfColors{nHist}, 'FaceColor', cvfColors{nHist}, 'LineWidth', 3);

   else
      % - Plot a stair series
      stairs(vfBins([1 1:end]), [0; cvfHist{nHist}(:); 0], '-', 'LineWidth', 3, 'Color', cvfColors{nHist});
   end
   
   hold on;
   
   % - Plot pips
   if (~isempty(cfPips) && ~isempty(cfPips{nHist}))
      plot(cfPips{nHist}, fPipLevel, 'v', ...
         'MarkerFaceColor', cvfColors{nHist}, 'MarkerEdgeColor', cvfColors{nHist}, 'LineWidth', 3, 'MarkerSize', 18);
   end
end

set(gca, 'FontSize', 32);
xlabel(strMeasureName);
ylabel('Count (%)');

% - Plot significance annotations
nPairIndex = 1;
for (nSeries1 = 1:numel(cvfData))
   for (nSeries2 = nSeries1+1:numel(cvfData))
      if (~isempty(cstrSig) && ~isempty(cstrSig{nPairIndex}) && ...
            ~isempty(cfPips) && ~isempty(cfPips{nSeries1}) && ~isempty(cfPips{nSeries2}))
         plot([cfPips{nSeries1} cfPips{nSeries2}], fPipLevel * 1.05 * [1 1], 'k-', 'LineWidth', 4);
         text(mean([cfPips{nSeries1} cfPips{nSeries2}]), fPipLevel * 1.075, cstrSig{nPairIndex}, 'HorizontalAlignment', 'center', 'FontSize', 48);
      end
      
      nPairIndex = nPairIndex + 1;
   end
end

set(gca, 'FontSize', 32, 'LineWidth', 3, 'TickLength', [0.01 0.01], 'TickDir', 'out', 'Box', 'off');

axis tight;
vfAxis = axis;
vfAxis([1 2]) = sign(vfBins([1 end])) .* abs(vfBins([1 end])) .* 1.05;
vfAxis([3 4]) = sign(vfAxis([3 4])) .* abs(vfAxis([3 4])) .* 1.05;
axis(vfAxis);

% - Add a title with annotation
set(gca, 'FontSize', 24);

title(['pips: ' sprintf('%g ', [cfPips{:}]) 'Ps: ' sprintf('%g ', vfP) 'Ns: ' sprintf('%d ', cellfun(@(c)(nnz(~isnan(c))), cvfData))]);
set(gca, 'FontSize', 32);

% - Restore hold state
if (~bHold)
   hold off;
end

% - Clear return argument, if not requested
if (nargout == 0)
   clear hFig;
end

% --- END of ComparisonHistN ---

%% -- Utility functions

% check_sig - FUNCTION Return the asterisk string to print for a given significance P values

function strSig = check_sig(fP)

if (fP < 0.001)
   strSig = '***';
elseif (fP < 0.01)
   strSig = '**';
elseif (fP < 0.05);
   strSig = '*';
else
   strSig = [];
end


% outputs - FUNCTION Select a particular output argument(s) from an anonymous function

function varargout = outputs(f, out_indices, varargin)

[outs{1:max(out_indices)}] = f(varargin{:});
varargout = outs(out_indices);

% --- END of ComparisonHistN.m ---
