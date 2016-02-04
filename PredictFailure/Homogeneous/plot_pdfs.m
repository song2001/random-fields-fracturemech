% Vincente Pericoli
% UC Davis


function plot_pdfs(samples, lstarIndex, lstars, saveFigures)
%PLOT_PDFS
% plot each sample's failure PDF for a given deterministic l*
%
% Inputs-
%   samples     = VGPy structure of samples, containing the CDF info
%   lstarIndex  = integer index identifying which l* is requested
%   lstars      = vector of lstars S.T. lstars(lstarIndex) = l*
%   saveFigures = optional boolean to request the figures be saved to a pdf
%                 file (default = False)
%

% set input defaults
if nargin < 4, saveFigures = false; end

% obtain the names & number of samples
sNames = fieldnames(samples);
nums = length(sNames);

% obtain material type
material = samples.(sNames{1}).material;

% set the results struct attribute name (based on lstarIndex)
str_res = sprintf('results%i',lstarIndex);

% loop through all samples
for s = 1:nums
    
    %
    % store data into new vars for convenience/readability
    %
    lkhoods    = samples.(sNames{s}).(str_res).lkhoods;
    failPDF    = samples.(sNames{s}).(str_res).failPDF;
    loadHist   = samples.(sNames{s}).(str_res).loadHist;
    failIndx   = samples.(sNames{s}).(str_res).failureIndex;
    distType   = samples.(sNames{s}).(str_res).distType;
    distParams = samples.(sNames{s}).(str_res).distParams;
    
    %
    % define title
    %
    
	% replace name underscores with spaces
    title_ = strrep(sNames{s},'_',' ');
    % provide info on distribution
    title_ = sprintf([title_, '\ndistribution =  %s'],distType);
    title_ = sprintf([title_, '\nmean =  %5.3f; stdev = %5.3f'], ...
                                    distParams(1),distParams(2));
    % provide the l*... switch format spec if l* is very small
	if lstars(lstarIndex) < 1e-4
        title_ = sprintf([title_,'\nl* = %4.2e'],lstars(lstarIndex));
    else
        title_ = sprintf([title_,'\nl* = %6.4f'],lstars(lstarIndex));
	end
    
    %
    % plot figure
    %
    H = figure;
	plot(loadHist,failPDF, loadHist,failPDF,'b*'); hold on;
    plot(loadHist(failIndx),lkhoods,'ro');
    % labels
    legend('PDF','Failure Observations', 'Location','Best');
    title(title_);
    xlabel('Load History');
    ylabel('Failure PDF');
    
    % if saveFigures is requested, save to a (temporary) .pdf file
    if saveFigures
        fileNames = strcat(sNames,'.pdf');
        saveas(H,fileNames{s},'pdf')
    end
end

%
% if savePDFs is requested, merge all temporary .pdf into a single file
%
if saveFigures
    % concatenate all files into a single .pdf
    filename = sprintf(['FailurePDFs_',material,'_','lstar%i.pdf'],lstarIndex);
    append_pdfs(filename,fileNames{:});
    % delete all temporary files
	delete(fileNames{:});
end

return
end