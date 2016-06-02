function tests = generateReportTest
  tests = functiontests(localfunctions);
end

function testGetAlgColors(testCase)
  % empty input should not generate error
  verifyEmpty(testCase, getAlgColors());
  % size of output should be n x 3
  for i = 0:2
    n = 10^i;
    verifySize(testCase, getAlgColors(1:n), [n, 3]);
  end
  % colors should be in range 0-255
  verifyGreaterThanOrEqual(testCase, getAlgColors(1:n), 0);
  verifyLessThanOrEqual(testCase, getAlgColors(1:n), 255);
end

function testCatEvalSet(testCase)
  % empty input should not generate error
  [evals, settings] = catEvalSet();
  verifyEmpty(testCase, evals);
  verifyEmpty(testCase, settings);
  % gain list of experiment folders containing scmaes_params.mat
  [~, expFolderList, ~] = expList();
  nExp = length(expFolderList);
  
  % test if there exists an experiment to report
  assumeNotEmpty(testCase, expFolderList, 'Cannot test generating report of one experiment. There is no folder containing experiment.');
  funcSet.BBfunc = [1 3 5 7 10];
  funcSet.dims = [2 3 20];
  [evals, settings] = catEvalSet(expFolderList(2), funcSet);
  % verify eval size
  [e1, e2, e3] = size(evals);
  verifyEqual(testCase, e1, length(funcSet.BBfunc));
  verifyEqual(testCase, e2, length(funcSet.dims));
  verifyEqual(testCase, e3, length(settings));
  
  % test if there exists more than 1 experiment to report
  assumeGreaterThanOrEqual(testCase, nExp, 2, 'Cannot test generating report of multiple experiment. There is only one folder containing experiment.');
  [evals, settings] = catEvalSet(expFolderList, funcSet);
  % verify eval size
  verifySize(testCase, evals, [length(funcSet.BBfunc), length(funcSet.dims), length(settings)]);
end

function testGenerateReport(testCase)
  % empty input should not generate error
  generateReport();
  % gain list of experiment folders containing scmaes_params.mat
  [folderList, expFolderList, expNameList] = expList();
  nExp = length(expFolderList);
  
  % test if there exists an experiment to report
  assumeNotEmpty(testCase, expFolderList, 'Cannot test generating report of one experiment. There is no folder containing experiment.');
  % verify existence of one-experiment report
  generateReport(expFolderList(1))
  verifyTrue(testCase, isdir(fullfile(expFolderList{1}, 'pproc')));
  verifyTrue(testCase, logical(exist(fullfile(expFolderList{1}, 'pproc', [expNameList{1}, '_report.m']), 'file')));
  
  % test if there exists more than 1 experiment to report
  assumeGreaterThanOrEqual(testCase, nExp, 2, 'Cannot test generating report of multiple experiment. There is only one folder containing experiment.');
  % verify existence of multi-experiment report
  generateReport(folderList)
  reportFile = cell(nExp, 1);
  for f = 1:nExp
    verifyTrue(testCase, isdir(fullfile(expFolderList{f}, 'pproc')));
    actualPP = dir(fullfile(expFolderList{f}, 'pproc'));
    actualPP = {actualPP(:).name};
    reportID = cellfun(@(x) ~isempty(strfind(x, [num2str(nExp), 'report'])), actualPP);
    % existence
    verifyTrue(testCase, any(reportID));
    reportFile{f} = actualPP(reportID);
  end
  uniqueReports = unique([reportFile{:}]);
  containReport = false(nExp, length(uniqueReports));
  for f = 1:nExp
    containReport(f,:) = cellfun(@(x) any(strcmp(x, reportFile{f})), uniqueReports);
  end
  % same report in all folders
  verifyTrue(testCase, any(all(containReport)));
  
  %TODO: publish generated report
  % problem - wrong paths due to unit test
  ppDir = fullfile(expFolderList{1}, 'pproc');
  addpath(ppDir)
  reportName = fullfile(ppDir, uniqueReports{find(all(containReport), 1, 'first')});
  publish(reportName)
end

function [folderList, expFolderList, expNameList] = expList()
  % gain list of experiment folders containing scmaes_params.mat
  actualFolder = pwd;
  cd(fullfile('..', '..'))
  scmaes_folder = pwd;
  cd(actualFolder)
  expFolder = fullfile(scmaes_folder, 'exp', 'experiments');
  folderNameList = dir(expFolder);
  folderNameList = {folderNameList([folderNameList(:).isdir]).name};
  folderNameList = folderNameList(3:end);
  folderList = cellfun(@(x) fullfile(expFolder, x), folderNameList, 'UniformOutput', false);
  expFolderID = logical(cellfun(@(x) exist(fullfile(x, 'scmaes_params.mat'), 'file'), folderList));
  expNameList = folderNameList(expFolderID);
  expFolderList = folderList(expFolderID);
end