import os
import json
import subprocess
from jobflow_remote import ConfigManager, JobController


def get_remote_run_dirs(function_name, flow_id_start, flow_id_end):
    cm = ConfigManager()
    jc = JobController.from_project(cm.get_project())
    run_dirs = []
    db_ids = [str(x) for x in range(flow_id_start, flow_id_end + 1)]
    while len(db_ids) > 0:
        db_id = db_ids.pop(0)
        flows_info_list = jc.get_flows_info(db_ids=db_id)
        if len(flows_info_list) == 0:
            # this id belongs to a deleted flow
            continue
        flows_info = flows_info_list[0]
        job_ids = flows_info.db_ids
        tmp = []
        for jid in job_ids:
            job_info = jc.get_jobs_info(db_ids=jid)[0]
            if jid in db_ids:
                idx = db_ids.index(jid)
                db_ids.pop(idx)
            if job_info.state.name != "COMPLETED":
                print(f"d_id={jid} is not finished, skipping!")
                continue
            if job_info.name == function_name:
                tmp.append(float(jid))
        if len(tmp) > 0:
            final_result_id = int(max(tmp))
            job_info = jc.get_jobs_info(db_ids=str(final_result_id))[0]
            run_dirs.append(job_info.run_dir)
    print(run_dirs)
    return run_dirs


def download_remote_jsons(
    file_name: str,
    file_dirs: list[str],
    cluster_name: str,
    store_dir: str | None = None
):
    if store_dir is None:
        store_dir = os.getcwd()
    for fdir in file_dirs:
        fdir = os.path.join(fdir, file_name)
        file = f"{cluster_name}:{fdir}"
        command = ['scp', file, store_dir]
        try:
            subprocess.run(command, check=True)
        except :
            print(f'can not download from {file}')
            continue
        tmp = os.path.join(store_dir, file_name)
        with open(tmp, 'r+') as fp:
            data = json.load(fp)
        elment = data['element']
        new_name = f'{elment}-{file_name}'
        new_name = os.path.join(store_dir, new_name)
        if os.path.exists(new_name):
            print(f'WARNING: {new_name} will be overwritten.')
        command = ['mv', tmp, new_name]
        print(f'{elment} is finished')
        subprocess.run(command, check=True)


if __name__ == "__main__":
    func_name = "export_result"
    file_name = 'eos_converge_results.json'
    store = '/home/wjing/PycharmProjects/pseudos_generation/ONCVPSP-PBE-SR-PDv0.6/delta1-v2'
    start = 276592
    end = 281157

    r_dir = get_remote_run_dirs(func_name, start, end)
    download_remote_jsons(file_name, r_dir, 'lucia', store_dir=store)
