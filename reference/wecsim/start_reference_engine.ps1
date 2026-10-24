$ErrorActionPreference = 'Stop'
$running = Get-Process MATLAB -ErrorAction SilentlyContinue
if ($running) { throw 'MATLAB already exists; reconcile before bootstrap' }
$startupCommand = "matlab.engine.shareEngine('CAPYTAINE_REFERENCE');set(0,'DefaultFigureVisible','off');"
$matlabProcess = Start-Process -FilePath matlab -ArgumentList '-nosplash','-nodesktop','-r',$startupCommand -WindowStyle Hidden -PassThru
[pscustomobject]@{launcher_pid=$matlabProcess.Id; launched_at=(Get-Date).ToUniversalTime().ToString('o'); host_name=$env:COMPUTERNAME; command=$startupCommand} | ConvertTo-Json | Set-Content -LiteralPath 'C:/Users/alex/topwater/artifacts/reference/reference-engine-bootstrap.json'
Get-Content -LiteralPath 'C:/Users/alex/topwater/artifacts/reference/reference-engine-bootstrap.json'
