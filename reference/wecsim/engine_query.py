"""Read-only shared MATLAB process/model inventory; no simulation or startup."""
import json
import matlab.engine

names = matlab.engine.find_matlab()
print(json.dumps({"shared_sessions": names}), flush=True)
if "CAPYTAINE_REFERENCE" not in names:
    raise SystemExit("Existing reference session is not visible")
engine = matlab.engine.connect_matlab("CAPYTAINE_REFERENCE")
try:
    print(engine.eval("jsonencode(struct('pid',feature('getpid'),'models',{find_system('type','block_diagram')},'variables',{who}))", nargout=1), flush=True)
    print(engine.eval("jsonencode(get_param(find_system('type','block_diagram'),'SimulationStatus'))", nargout=1), flush=True)
finally:
    engine.exit()
