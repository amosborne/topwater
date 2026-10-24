function fix_imported_dispatch()
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
if bdIsLoaded('minimum_reference')
    assert(strcmp(get_param('minimum_reference','SimulationStatus'),'stopped'));
    close_system('minimum_reference',0);
end
lib = 'WECSim_Lib_Body_Elements';
load_system(lib);
set_param(lib,'Lock','off');
parent = [lib '/Rigid Body/Hydrodynamic Body/Wave Diffraction and Excitation Force Calculation'];
mask = Simulink.Mask.get(parent);
old = mask.Initialization;
line = 'nWaves = size(waveGroup,2);';
assert(numel(strfind(old,line)) == 1 && ~contains(old,'importedExcitationBranches'));
% Edit this authored mask initialization through its native API; the original
% branch-creation loop remains verbatim for all non-imported wave modes.
mask.Initialization = strrep(old,line,'nWaves = importedExcitationBranches(cb,waveGroup);');
assert(contains(mask.Initialization,'nWaves = importedExcitationBranches(cb,waveGroup);'));
set_param(lib,'Lock','on');
save_system(lib);
close_system(lib,0);
end
