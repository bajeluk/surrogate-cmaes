
expdir = 'saACMESlambdaRevMinIter3v2';
nexp = 60;
if (1)
    
    
    id1 = 1;
    
    for iexp=1:nexp

        file = fopen(['script' num2str(iexp)],'w');
        dir = ['/users/tao/ilya/MatLab/ACMpaper2/' expdir];
     %   fprintf(file, ['rm -r ' dir '/res' num2str(iexp) '.txt\n']);
     %   fprintf(file, ['nohup /special/matlab/bin/matlab -nojvm -nodisplay <' dir '/exp' num2str(iexp) '.m> ' dir '/res' num2str(iexp) '.txt &']);
        fprintf(file, ['nohup octave -q ' dir '/exp' num2str(iexp) '.m ' dir '/res' num2str(iexp) '.txt &']);
        fclose(file);

    end;

end;
if (1)
    
    ntipi = 8;
    id1 = 1;
    iexp = 1;
    for itipi=1:ntipi
        file = fopen(['tipi' num2str(itipi)],'w');
        for irunOnTipi=1:8
            dir = ['/users/tao/ilya/MatLab/ACMpaper2/' expdir];
         %   fprintf(file, ['rm -r ' dir '/res' num2str(iexp) '.txt\n']);
         %   fprintf(file, ['nohup /special/matlab/bin/matlab -nojvm -nodisplay <' dir '/exp' num2str(iexp) '.m> ' dir '/res' num2str(iexp) '.txt &']);
            fprintf(file, ['nohup octave -q ' dir '/exp' num2str(iexp) '.m ' dir '/res' num2str(iexp) '.txt &\n']);
            iexp = iexp + 1;
        end;
        fclose(file);
    end;

end;
z = 0;