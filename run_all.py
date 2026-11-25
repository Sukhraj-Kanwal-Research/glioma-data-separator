import subprocess
import sys

scripts = [
    r"C:\Users\Owner\Desktop\prep_qc_and_log1p.py",
    r"C:\Users\Owner\Desktop\run_pca.py",
    r"C:\Users\Owner\Desktop\run_umap_cluster.py",
    r"C:\Users\Owner\Desktop\find_markers.py",
]

python_exe = sys.executable

for script in scripts:
    print(f"\n=== Running: {script} ===")
    result = subprocess.run([python_exe, script], check=False)
    if result.returncode != 0:
        print(f"[!] Stopped because {script} exited with code {result.returncode}")
        sys.exit(result.returncode)

print("\n[SAVED] Full pipeline completed cleanly!")