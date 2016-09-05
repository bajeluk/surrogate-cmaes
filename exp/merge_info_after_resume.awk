BEGIN   { RS="funcId"; ORS=""
          algId = ""; }
        { if (length($0) > 0) {
            if (algId != $11) {
              # algId has changed
              if (length(algId) > 0)
                printf "\n";    # technical newline because ORS == ""
              algId = $11;
              file = $13;
              sub(",$", "", file);
              sub("\n$", "");   # strip the last newline -- we can continue on the same line

              if (match(algId, "_cmaes'") && ! match($13, "exp-01_f"))
              {
                newCmaesFile = file;
                gsub("exp-0[0-9]_f", "exp-01_f", newCmaesFile);
                gsub("exp-0[0-9]_f", "exp-01_f", $0);
                cmd = "mv " file " " newCmaesFile;
                system(cmd);
                gsub("\\.dat", ".tdat", cmd);
                system(cmd);
                gsub("\\.tdat", ".rdat", cmd);
                system(cmd);
              }
              printf "funcId%s", $0;
            }
            else {
              # algId stays the same...
              # copy the contents of files into the first file (when algId changed)
              newFile = $13;
              sub(",$", "", newFile);
              cmd = "cat " newFile " >> " file;
              system(cmd);
              system("rm " newFile);
              gsub("\\.dat", ".tdat", cmd);
              system(cmd);
              gsub("\\.dat", ".tdat", newFile);
              system("rm " newFile);
              gsub("\\.tdat", ".rdat", cmd);
              system(cmd);
              gsub("\\.tdat", ".rdat", newFile);
              system("rm " newFile);
              # print "CMD: " cmd
              gsub("^.*data_f[^,]*", "");
              sub("\n$", "")
              printf "%s", $0;
              # printf "%s", $0;
            }
          }
        }
END     { printf "\n"; }
