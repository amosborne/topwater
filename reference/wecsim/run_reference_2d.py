"""Telemetry augmentation then the public shared rendering CLI; no plotting here."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from augment_scene_telemetry import augment,read,sha

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    c=read(args.config);started=time.monotonic();qa=ROOT/c['qa_directory'];qa.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,UV_PROJECT_ENVIRONMENT=str(Path(sys.executable).resolve().parents[1]),VIRTUAL_ENV=str(Path(sys.executable).resolve().parents[1]),UV_CACHE_DIR=str(ROOT/'artifacts/tools/uv-cache'),PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    subprocess.run([sys.executable,'-m','pytest','-q','--basetemp',str(qa/'pytest-temp'),'tests/test_visualization.py','tests/test_package.py'],cwd=ROOT,env=env,check=True)
    proof=augment(c);(qa/'augmentation.json').write_text(json.dumps(proof,indent=2)+'\n')
    command=[shutil.which('uv'),'run','--no-sync','topwater','render',str(ROOT/c['augmented_scene']),str(ROOT/c['output'])]
    process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    print(dict(public_cli=command,cli_pid=process.pid,parent_pid=os.getpid()),flush=True)
    stdout,stderr=process.communicate();(qa/'stdout.txt').write_text(stdout);(qa/'stderr.txt').write_text(stderr)
    assert process.returncode==0,f'Public CLI failed {process.returncode}: {stderr}'
    report=json.loads(stdout);(qa/'render-report.json').write_text(json.dumps(report,indent=2)+'\n')
    movie=ROOT/c['output']
    probe=subprocess.run([shutil.which('ffprobe'),'-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration','-of','json',str(movie)],check=True,capture_output=True,text=True)
    metadata=json.loads(probe.stdout);stream=metadata['streams'][0]
    assert int(stream['nb_read_frames'])==c['frames'] and stream['width']==800 and stream['height']==800
    assert stream['r_frame_rate']=='25/2' and abs(float(stream['duration'])-c['duration_s'])<1e-5
    (qa/'ffprobe.json').write_text(json.dumps(metadata,indent=2)+'\n')
    subprocess.run([shutil.which('ffmpeg'),'-v','error','-i',str(movie),'-vf',r'select=eq(n\,18)','-frames:v','1',str(qa/'poster.png')],check=True)
    assert all(sha(ROOT/p)==digest for p,digest in c['input_sha256'].items())
    (qa/'verification.json').write_text(json.dumps(dict(complete=True,elapsed_s=time.monotonic()-started,mp4_sha256=sha(movie),metadata=metadata,augmentation=proof,render_report=report),indent=2)+'\n')

if __name__=='__main__':main()
