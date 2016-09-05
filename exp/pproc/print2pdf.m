function print2pdf(handle, pdfname, overwrite)
% Print figures to pdf in format displayed on the screen.
%
% Input:
%   handle    - handles of figures | cell-array
%   pdfname   - string or cell array of strings containing resulting 
%               filenames
%   overwrite - automatically overwrite without asking | boolean | double

  if nargin < 1
    help print2pdf
    return
  end

  nFig = length(handle);
  if nargin < 2
    for f = 1:nFig
      pdfname{f} = ['figure', num2str(f), '.pdf'];
    end
  end
  % pass handle input to cell-array
  if ishandle(handle)
    handle = num2cell(handle);
  end
  if ischar(pdfname)
    pdfname = {pdfname};
  end

  nNames = length(pdfname);
  % check if names end with .pdf
  for f = 1:nNames
    if ~strcmp(pdfname{f}(end-3:end),'.pdf')
      pdfname{f} = [pdfname{f},'.pdf'];
    end
  end

  % check if there is enough names
  if nNames ~= nFig
    if nNames == 1
      fprintf('Creating names:\n\n')
      pdfname(1:nFig) = pdfname;
      for f = 1:nFig
        pdfname{f} = [pdfname{f}(1:end-4),num2str(f),pdfname{f}(end-3:end)];
        fprintf('%s\n',pdfname{f})
      end
    else
    error('Numbers of figures and names does not agree!')
    end
  end

  % count existing files
  existingPDFs = [];
  for f = 1:nFig
    pdfFolderId = strfind(pdfname{f}, filesep);
    pdfFolder = pdfname{f}(1:pdfFolderId(end) - 1);
    if ~exist(pdfFolder, 'dir')
      mkdir(pdfFolder)
    elseif exist(pdfname{f},'file')
      existingPDFs(end+1) = f;            
    end
  end
  NexistPDF = length(existingPDFs);

  if nargin < 3
    % ask for overwritting previous images (if there were any)
    if any(existingPDFs)
        answer = questdlg(['Overwrite pdfs for ',num2str(NexistPDF),' files?'],'Overwritting old files','Overwrite','No','Overwrite');
        if strcmp('Overwrite',answer)
            overwrite = 1;
        else
            overwrite = 0;
        end
    else
        overwrite = 1;
    end
  end 

  % print plot to pdf
  if overwrite
    for f = 1:nFig
      set(handle{f},'PaperPositionMode','auto')
      print(handle{f},'-dpdf','-r0',pdfname{f});
    end
  else
    for f = 1:NexistPDF
      set(handle{existingPDFs(f)},'PaperPositionMode','auto')
      print('-dpdf','-r0',pdfname{existingPDFs(f)});
    end
  end

end