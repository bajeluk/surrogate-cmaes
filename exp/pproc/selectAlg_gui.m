function setID = selectAlg_gui(expSettings)
% selectAlg_gui Select an algorithm settings from the listbox, then
% add it to list to report. 

  % initialization
  nSettings = length(expSettings);
  setID = 1 : nSettings;
  algNames = arrayfun(@(x) ['ALG', num2str(x)], 1:nSettings, 'UniformOutput', false);
  % find differences
  [dFields, dValues] = difField(expSettings);
  
  % component sizes
  marg  = 20;
  margIn = 60;
  listBoxSize = [100, 150];
  buttonSize  = [70,  25];
  setSize    = [2*listBoxSize(1) + 2*margIn + buttonSize(1), 105];
  formSize    = [2*marg + setSize(1), 4*marg + listBoxSize(2) + setSize(2) + buttonSize(2)];
  
  % component positions
  allPosX = marg;
  allPosY = 3*marg + buttonSize(2) + setSize(2);
  chosenPosX = marg +2*margIn + listBoxSize(1) + buttonSize(1);
  chosenPosY = 3*marg + setSize(2) + buttonSize(2);
  addPosX = marg + margIn + listBoxSize(1);
  addPosY = 3*marg + setSize(2) + buttonSize(2) + (3/4 + 1/8)*listBoxSize(2) - buttonSize(2)/2;
  remPosX = addPosX;
  remPosY = addPosY - 1/4*listBoxSize(2);
  addAllPosX = addPosX;
  addAllPosY = remPosY - 1/4*listBoxSize(2);
  remAllPosX = addPosX;
  remAllPosY = addAllPosY - 1/4*listBoxSize(2);
  setPosX = marg;
  setPosY = 2*marg + buttonSize(2);
  okPosX = 1/4*formSize(1) - buttonSize(1)/2;
  okPosY = marg;
  cancelPosX = 3/4*formSize(1) - buttonSize(1)/2;
  cancelPosY = okPosY;

  % Create and then hide the UI as it is being constructed.
  handle.f = figure('Visible', 'off', 'Position', [360, 500, formSize]);

  % Construct the components.
  handle.AllSettings  = uicontrol('Style', 'listbox', ...
                   'Position', [ allPosX, allPosY, listBoxSize], ...
                   'Callback', @allSettings_Callback, ...
                   'String', algNames);
  handle.ChosenSettings  = uicontrol('Style', 'listbox', ...
                      'Position', [ chosenPosX, chosenPosY, listBoxSize], ...
                      'Callback', @chosenSettings_Callback, ...
                      'String', []);
  handle.Add       = uicontrol('Style', 'pushbutton', ...
                'Position', [addPosX, addPosY, buttonSize], ...
                'Callback', @addButton_Callback, ...
                'String', '>');
  handle.Remove    = uicontrol('Style', 'pushbutton', ...
                'Position', [remPosX, remPosY, buttonSize], ...
                'Callback', @removeButton_Callback, ...
                'String', '<');
  handle.AddAll    = uicontrol('Style', 'pushbutton', ...
                'Position', [addAllPosX, addAllPosY, buttonSize], ...
                'Callback', @addAllButton_Callback, ...
                'String', '>>');
  handle.RemoveAll = uicontrol('Style', 'pushbutton', ...
                'Position', [remAllPosX, remAllPosY, buttonSize], ...
                'Callback', @removeAllButton_Callback, ...
                'String', '<<');
  handle.OK        = uicontrol('Style', 'pushbutton', ...
                'Position', [okPosX, okPosY, buttonSize], ...
                'Callback', {@okButton_Callback, handle}, ...
                'String', 'OK');
  handle.Cancel    = uicontrol('Style', 'pushbutton', ...
                'Position', [cancelPosX, cancelPosY, buttonSize], ...
                'Callback', @cancelButton_Callback, ...
                'String', 'Cancel');
  handle.Settings  = uicontrol('Style', 'edit', ...
                'Position', [setPosX, setPosY, setSize], ...
                'Callback', @settingsButton_Callback, ...
                'String', returnDiffField(dFields, dValues, 1), ...
                'Max', 30, ...
                'FontName', 'Monospaced', ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'left');

% assign a name to appear in the window title
handle.f.Name = 'Choose algorithms to report';
% move the window to the center of the screen
movegui(handle.f, 'center')
% make the window visible
handle.f.Visible = 'on';

end

function s = returnDiffField(dField, dValue, setNumber)
% return string containing different fields in settings
  s = cell(length(dField), 1);
  for f = 1:length(dField)
    s{f} = sprintf('%s = %s;', dField{f}, printStructure(dValue{f, setNumber}, 1, 'Format', 'value'));
  end
end

function okButton_Callback(source, eventdata, handle)
  close
end

function cancelButton_Callback(source, eventdata)
  close
end