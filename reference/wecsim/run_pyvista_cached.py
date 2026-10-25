"""Cached standard-CLI rerender only; no wave preparation or simulation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    c=read(args.config);started=time.monotonic();assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    qa=ROOT/c['qa_directory'];qa.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,UV_PROJECT_ENVIRONMENT=str(Path(sys.executable).resolve().parents[1]),VIRTUAL_ENV=str(Path(sys.executable).resolve().parents[1]),UV_CACHE_DIR=str(ROOT/'artifacts/tools/uv-cache'),PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VTK_SMP_MAX_THREADS='1')
    command=[shutil.which('uv'),'run','--no-sync','topwater','render',str(ROOT/c['scene']),str(ROOT/c['output']),'--settings',str(ROOT/c['settings'])]
    process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    print(dict(public_cli=command,cli_pid=process.pid,parent_pid=os.getpid()),flush=True)
    stdout,stderr=process.communicate();(qa/'stdout.txt').write_text(stdout);(qa/'stderr.txt').write_text(stderr)
    assert process.returncode==0,f'Public renderer CLI failed {process.returncode}: {stderr}'
    report=json.loads(stdout);write(qa/'render-report.json',report)
    movie=ROOT/c['output']
    result=subprocess.run([shutil.which('ffprobe'),'-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration','-of','json',str(movie)],capture_output=True,text=True,check=True)
    metadata=json.loads(result.stdout);stream=metadata['streams'][0]
    assert int(stream['nb_read_frames'])==c['frames'] and stream['width']==800 and stream['height']==800
    assert stream['r_frame_rate']=='25/2' and abs(float(stream['duration'])-c['duration_s'])<1e-5
    write(qa/'ffprobe.json',metadata)
    subprocess.run([shutil.which('ffmpeg'),'-v','error','-i',str(movie),'-vf',r'select=eq(n\,18)','-frames:v','1',str(qa/'poster.png')],check=True)
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    write(qa/'verification.json',dict(complete=True,elapsed_s=time.monotonic()-started,mp4_sha256=sha(movie),scene_sha256=sha(ROOT/c['scene']),metadata=metadata,render_report=report,scope='Cached-only standardPyVista rerender; no field/physics/preparation/test/install changes'))

if __name__=='__main__':main()
