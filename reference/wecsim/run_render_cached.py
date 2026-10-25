"""Cached-only public rendering CLI plus explicit reference timing and QA."""
import argparse,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')

def main():
    started=time.perf_counter();parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args();c=read(args.config)
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    qa=ROOT/c['qa_directory'];qa.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,UV_PROJECT_ENVIRONMENT=str(Path(sys.executable).resolve().parents[1]),VIRTUAL_ENV=str(Path(sys.executable).resolve().parents[1]),UV_CACHE_DIR=str(ROOT/'artifacts/tools/uv-cache'),PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    focused_tests_seconds=0.
    if c.get('run_focused_tests',True):
        test_started=time.perf_counter()
        subprocess.run([sys.executable,'-m','pytest','-q','--basetemp',str(qa/'pytest-temp'),'tests/test_visualization.py::test_header_columns_stay_fixed_across_signs_and_digit_changes','tests/test_visualization.py::test_opaque_footprint_uses_projected_input_geometry'],cwd=ROOT,env=env,check=True)
        focused_tests_seconds=time.perf_counter()-test_started
    command=[shutil.which('uv'),'run','--no-sync','topwater','render',str(ROOT/c['scene']),str(ROOT/c['output'])]
    cli_started=time.perf_counter();process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    print(dict(public_cli=command,cli_pid=process.pid,parent_pid=os.getpid()),flush=True)
    stdout,stderr=process.communicate();cli_wall=time.perf_counter()-cli_started
    (qa/'stdout.txt').write_text(stdout);(qa/'stderr.txt').write_text(stderr)
    assert process.returncode==0,f'Public CLI failed {process.returncode}: {stderr}'
    report=json.loads(stdout);write(qa/'render-report.json',report)
    qa_started=time.perf_counter();movie=ROOT/c['output']
    probe=subprocess.run([shutil.which('ffprobe'),'-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration','-of','json',str(movie)],capture_output=True,text=True,check=True)
    metadata=json.loads(probe.stdout);s=metadata['streams'][0]
    assert int(s['nb_read_frames'])==c['frames'] and s['width']==800 and s['height']==800
    assert s['r_frame_rate']=='25/2' and abs(float(s['duration'])-c['duration_s'])<1e-5
    write(qa/'ffprobe.json',metadata)
    subprocess.run([shutil.which('ffmpeg'),'-v','error','-i',str(movie),'-vf',r'select=eq(n\,18)','-frames:v','1',str(qa/'poster.png')],check=True)
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    write(qa/'verification.json',dict(complete=True,runner_wall_seconds=time.perf_counter()-started,external_cli_wall_seconds=cli_wall,focused_tests_seconds=focused_tests_seconds,qa_wall_seconds=time.perf_counter()-qa_started,mp4_sha256=sha(movie),scene_sha256=sha(ROOT/c['scene']),metadata=metadata,render_report=report,scope='Exact already-telemetryaugmentedNPZ cached-only; no augmentation/waveprep/physics'))

if __name__=='__main__':main()
