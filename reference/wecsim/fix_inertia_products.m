function fix_inertia_products()
root = 'C:/Users/alex/topwater';
addpath(genpath(fullfile(root,'artifacts/dependencies/topwater-wec-sim/source')));
lib = 'WECSim_Lib_Body_Elements';
load_system(lib);
set_param(lib,'Lock','off');
blocks = {'Flex Body/Structure/Body Properties', ...
    'Rigid Body/Drag Body/Structure/Body Properties', ...
    'Rigid Body/Hydrodynamic Body/Structure/Body Properties'};
for i = 1:3
    block = [lib '/' blocks{i}];
    old = get_param(block,'ProductsOfInertia');
    assert(~contains(old,'[3 2 1]'),'Library already reordered');
    assert(contains(old,'inertiaProducts') || contains(old,'adjustedInertiaProducts'));
    % public/internal WEC-Sim [Ixy Ixz Iyz] to File Solid [Iyz Izx Ixy]
    updated = [old '([3 2 1])'];
    set_param(block,'ProductsOfInertia',updated);
    assert(strcmp(get_param(block,'ProductsOfInertia'),updated));
end
set_param(lib,'Lock','on');
save_system(lib);
close_system(lib,0);
check_inertia_mapping();
end
