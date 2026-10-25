classdef WorkflowEngineGuard < handle
    % Restore shared-session settings and close only this run's unique model.
    properties (Access = private)
        state
        model
        restored = true
    end
    methods
        function obj = WorkflowEngineGuard(model,originalPath)
            assert(~bdIsLoaded(model),'Workflow model name is already loaded');
            obj.model = model;
            cfg = Simulink.fileGenControl('getConfig');
            % Copy values: getConfig returns a HANDLE which later setters mutate.
            obj.state = struct('directory',pwd,'path',originalPath, ...
                'diary',get(0,'Diary'),'diaryFile',get(0,'DiaryFile'), ...
                'warnings',warning('query','all'),'threads',maxNumCompThreads, ...
                'cacheFolder',cfg.CacheFolder,'codeGenFolder',cfg.CodeGenFolder, ...
                'codeGenFolderStructure',cfg.CodeGenFolderStructure, ...
                'figureVisible',get(groot,'DefaultFigureVisible'));
            obj.restored = false;
        end

        function prepare(obj)
            % Preserve the user's desktop/minimized state; never activate it.
            diary off
            set(groot,'DefaultFigureVisible','off');
        end

        function report = restore(obj)
            if obj.restored
                report = struct('restored',true,'failures',{{}});
                return
            end
            s = obj.state;
            failures = {};
            actions = {@() close_model(obj.model), @() diary('off'), ...
                @() cd(s.directory), ...
                @() Simulink.fileGenControl('set','CacheFolder',s.cacheFolder, ...
                    'CodeGenFolder',s.codeGenFolder, ...
                    'CodeGenFolderStructure',s.codeGenFolderStructure, ...
                    'keepPreviousPath',false,'createDir',false), ...
                @() path(s.path), @() maxNumCompThreads(s.threads), ...
                @() set(groot,'DefaultFigureVisible',s.figureVisible), ...
                @() warning(s.warnings), @() restore_diary(s)};
            labels = {'owned model','workflow diary','working directory', ...
                'file generation folders','search path','thread limit', ...
                'figure visibility default','warnings','original diary'};
            for k = 1:numel(actions)
                try
                    actions{k}();
                catch failure
                    failures{end+1} = [labels{k} ': ' failure.message];
                end
            end
            obj.restored = isempty(failures);
            report = struct('restored',obj.restored,'failures',{failures});
        end

        function delete(obj)
            % Fallback for errors/interrupts before Python's explicit finally.
            if ~obj.restored
                report = obj.restore();
                if ~report.restored
                    warning('Topwater shared-engine restoration incomplete: %s', ...
                        strjoin(report.failures,'; '));
                end
            end
        end
    end
end

function close_model(model)
if bdIsLoaded(model)
    if ~strcmp(get_param(model,'SimulationStatus'),'stopped')
        set_param(model,'SimulationCommand','stop');
    end
    close_system(model,0);
end
end

function restore_diary(state)
set(0,'DiaryFile',state.diaryFile);
diary(state.diary);
end
