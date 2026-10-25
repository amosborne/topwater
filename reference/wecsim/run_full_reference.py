"""Full accepted-cache/border data assembly then the unchanged public 2D CLI."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from prepare_full_reference_scene import prepare_full_scene,read,sha,write

def main():
    started=time.perf_counter();parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);a=parser.parse_args();c=read(a.config)
    qa=ROOT/c['qa_directory'];qa.mkdir(parents=True,exist_ok=False)
    proof=prepare_full_scene(c);write(qa/'preparation.json',proof);print(json.dumps(proof),flush=True)
    env=dict(os.environ,UV_PROJECT_ENVIRONMENT=str(Path(sys.executable).resolve().parents[1]),VIRTUAL_ENV=str(Path(sys.executable).resolve().parents[1]),UV_CACHE_DIR=str(ROOT/'artifacts/tools/uv-cache'),PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    command=[shutil.which('uv'),'run','--no-sync','topwater','render',str(ROOT/c['scene_output']),str(ROOT/c['output'])]
    before=time.perf_counter();process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    print(dict(public_cli=command,cli_pid=process.pid,parent_pid=os.getpid()),flush=True)
    stdout,stderr=process.communicate();cli_wall=time.perf_counter()-before
    (qa/'stdout.txt').write_text(stdout);(qa/'stderr.txt').write_text(stderr)
    assert process.returncode==0,f'Public CLI failed {process.returncode}: {stderr}'
    report=json.loads(stdout);write(qa/'render-report.json',report)
    before=time.perf_counter();movie=ROOT/c['output']
    probe=subprocess.run([shutil.which('ffprobe'),'-v','error','-count_frames','-select_streams','v:0','-show_entries','stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration','-of','json',str(movie)],check=True,capture_output=True,text=True)
    metadata=json.loads(probe.stdout);s=metadata['streams'][0]
    assert int(s['nb_read_frames'])==c['frames'] and s['width']==800 and s['height']==800
    assert s['r_frame_rate']=='25/2' and abs(float(s['duration'])-c['duration_s'])<1e-5
    write(qa/'ffprobe.json',metadata)
    for index in c['qa_frames']:
        subprocess.run([shutil.which('ffmpeg'),'-v','error','-i',str(movie),'-vf',f'select=eq(n\\,{index})','-frames:v','1',str(qa/f'frame-{index:04d}.png')],check=True)
    assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
    assert sha(ROOT/c['scene_output'])==proof['scene_sha256']
    write(qa/'verification.json',dict(complete=True,runner_wall_seconds=time.perf_counter()-started,preparation_seconds=proof['preparation_seconds'],external_cli_wall_seconds=cli_wall,qa_wall_seconds=time.perf_counter()-before,mp4_sha256=sha(movie),scene_sha256=proof['scene_sha256'],metadata=metadata,preparation=proof,render_report=report))

if __name__=='__main__':main()
