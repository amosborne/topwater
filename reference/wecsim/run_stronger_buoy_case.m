%%
%[text] Stronger observed-buoy case: data/path entry only, shared native physics.
%%
strongerPreviousRoot = getenv('TOPWATER_BUOY_ROOT');
strongerPreviousTrain = getenv('TOPWATER_BUOY_TRAIN');
setenv('TOPWATER_BUOY_ROOT','C:/Users/alex/topwater/artifacts/reference/stronger-buoy-20261024');
setenv('TOPWATER_BUOY_TRAIN','train-attempt1');
try
    run('C:/Users/alex/topwater/reference/wecsim/run_calm_buoy_case.m');
catch caseError
    setenv('TOPWATER_BUOY_ROOT',strongerPreviousRoot);
    setenv('TOPWATER_BUOY_TRAIN',strongerPreviousTrain);
    rethrow(caseError)
end
setenv('TOPWATER_BUOY_ROOT',strongerPreviousRoot);
setenv('TOPWATER_BUOY_TRAIN',strongerPreviousTrain);
%[appendix]{"version":"1.0"}
