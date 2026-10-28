# Start the shared WEC-Sim engine with temporary files outside the checkout.
$ErrorActionPreference = 'Stop'
$running = Get-Process MATLAB -ErrorAction SilentlyContinue
if ($running) { throw 'MATLAB already exists; reconcile before starting another engine' }

$profileDirectory = [Environment]::GetFolderPath('UserProfile')
$runtimeDirectory = Join-Path $profileDirectory ('.cache/topwater/matlab/' + [guid]::NewGuid().ToString('N'))
$temporaryDirectory = Join-Path $runtimeDirectory 'tmp'
$cacheDirectory = Join-Path $runtimeDirectory 'cache'
$codeGenerationDirectory = Join-Path $runtimeDirectory 'codegen'
New-Item -ItemType Directory -Path $temporaryDirectory, $cacheDirectory, $codeGenerationDirectory -Force | Out-Null

# MATLAB inherits these values. Restore this PowerShell process after launch.
$originalTemporary = [Environment]::GetEnvironmentVariable('TMP', 'Process')
$originalTemp = [Environment]::GetEnvironmentVariable('TEMP', 'Process')
$matlabCache = $cacheDirectory.Replace("'", "''")
$matlabCodeGeneration = $codeGenerationDirectory.Replace("'", "''")
$startupCommand = "set(groot,'DefaultFigureVisible','off');" +
    "Simulink.fileGenControl('set','CacheFolder','$matlabCache','CodeGenFolder','$matlabCodeGeneration','createDir',true);" +
    "matlab.engine.shareEngine('CAPYTAINE_REFERENCE');"
try {
    [Environment]::SetEnvironmentVariable('TMP', $temporaryDirectory, 'Process')
    [Environment]::SetEnvironmentVariable('TEMP', $temporaryDirectory, 'Process')
    $launchArguments = @('-nosplash', '-nodesktop', '-sd', ('"{0}"' -f $runtimeDirectory),
        '-r', ('"{0}"' -f $startupCommand))
    $matlabProcess = Start-Process -FilePath matlab -ArgumentList $launchArguments -WindowStyle Hidden -PassThru
} finally {
    [Environment]::SetEnvironmentVariable('TMP', $originalTemporary, 'Process')
    [Environment]::SetEnvironmentVariable('TEMP', $originalTemp, 'Process')
}

# Sharing completes during startup; the adapter checks engine visibility.
[pscustomobject]@{
    launcher_pid = $matlabProcess.Id
    engine = 'CAPYTAINE_REFERENCE'
    runtime_directory = $runtimeDirectory
}
