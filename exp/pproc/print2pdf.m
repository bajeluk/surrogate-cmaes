function print2pdf(handle, pdfname, overwrite)
% Print figures to pdf in format displayed on the screen.
%
% Input:
%   handle    - array of handles of figures
%   pdfname   - string or cell array of strings containing resulting 
%               filenames
%   overwrite - automatically overwrite without asking | boolean | double

  if nargin == 0
    help print2pdf
    return
  end

  Nfig = length(handle);
  if ischar(pdfname)
    pdfname = {pdfname};
  end

  Nnames = length(pdfname);
  % check if names end with .pdf
  for f = 1:Nnames
    if ~strcmp(pdfname{f}(end-3:end),'.pdf')
      pdfname{f} = [pdfname{f},'.pdf'];
    end
  end

  % check if there is enough names
  if Nnames ~= Nfig
    if Nnames == 1
      fprintf('Creating names:\n\n')
      pdfname(1:Nfig) = pdfname;
      for f = 1:Nfig
        pdfname{f} = [pdfname{f}(1:end-4),num2str(f),pdfname{f}(end-3:end)];
        fprintf('%s\n',pdfname{f})
      end
    else
    error('Numbers of figures and names does not agree!')
    end
  end

  % count existing files
  existingPDFs = [];
  for f = 1:Nfig
    if exist(pdfname{f},'file')
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
    for f = 1:Nfig
      set(handle(f),'PaperPositionMode','auto')
      print(handle(f),'-dpdf','-r0',pdfname{f});
    end
  else
    for f = 1:NexistPDF
      set(handle(existingPDFs(f)),'PaperPositionMode','auto')
      print('-dpdf','-r0',pdfname{existingPDFs(f)});
    end
  end

end