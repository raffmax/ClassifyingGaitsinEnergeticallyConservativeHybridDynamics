function res = deleteAppData(identifier)
    %DELETEAPPDATA deletes AppData cache or the folders called "identifier" in AppData
    %   Function will not ask for user input or do anything if no AppData folder exists.
    %   Input: none - delete complete AppData cache
    %       OR identifier (string) of an Mechanical Model (corresponds to
    %          the subfolder name in AppData)
    %   Output: +2 if AppData folder gets deleted
    %           +1 if subfolder gets deleted
    %           -1 if nothing gets deleted
    %           -2 if AppData or AppData/identifier doesn't exist
    
    
    switch nargin
        case 0 % no subfolder name given -> delete complete AppData
            
            ad_pth = what('AppData'); % lists all appdata/model paths
            if ~isempty(ad_pth)
                disp('---- Folder[s] listed below ----')
                for idx=1:length(ad_pth)
                    disp(ad_pth(idx).path);
                end
                disp('---- Folder[s] listed above ----')
                validInput = false;
                while ~validInput 
                    % ask for input
                    del = input('Delete AppData? [y/n] \n', 's');
                    if any(strcmpi(del, {'y', 'n', 'all', 'yes', 'no'}))
                        validInput = true;
                    else
                        disp('Invalid input. Allowed input: [y/n]');
                    end
                    
                    % - delete - 
                    if any(strcmpi(del, {'y', 'yes'})) 
                        addpath(genpath(ad_pth.path))
                        rmpath(genpath(ad_pth.path))
                        rmdir(ad_pth(idx).path, 's');
                        res = 2;
                    else
                        res = -1; 
                    end
                end
            else
                res = -2; % AppData folder doesn't exist
            end

        case 1 % ask if subfolder or complete AppData should be cleard
            ad_dir = ['AppData', filesep, identifier];
            ad_pth = what(ad_dir); % lists all appdata/model paths
            if ~isempty(ad_pth)
                disp('---- Folder[s] listed below ----')
                for idx=1:length(ad_pth)
                    disp(ad_pth(idx).path);
                end
                 disp('---- Folder[s] listed above ----')
                validInput = false;
                while ~validInput
                    del = input('Delete above mentioned folder[s] from AppData? [y/n/all] \nNote: \"all\" deletes all folders in the AppData cache. \n', 's');
                    % check for valid input
                    if any(strcmpi(del, {'y', 'n', 'all', 'yes', 'no'})) 
                        validInput = true;
                    else 
                        disp('Invalid input. Allowed input: [y/n/all]');
                    end
                    
                    % delete
                    if strcmpi(del, 'y') % del==y==Y
                        addpath(ad_pth(idx).path)
                        rmpath(ad_pth(idx).path)
                        for idx=1:length(ad_pth)
                            rmdir(ad_pth(idx).path, 's');
                        end
                        res = +1;
                    elseif strcmpi(del, 'all') % delete complete AppData folder[s]
                        disp('Deleting AppData folder[s] ...');
                        ad_pth = what('AppData');
                        addpath(genpath(ad_pth.path))
                        rmpath(genpath(ad_pth.path))
                        rmdir(ad_pth(idx).path, 's');
                        res = +2;
                    else 
                        res = -1;
                    end
                    
                end
            else
                res = -2; % AppData or AppData/identifier doesn't exist
            end        
    end
end

