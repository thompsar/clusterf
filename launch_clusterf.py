import os
import sys
import subprocess

# Set PYTHONPATH to current directory (parent of clusterf/)
project_root = os.path.abspath(os.path.dirname(__file__))
env = os.environ.copy()
env["PYTHONPATH"] = project_root

# Call `panel serve main.py --show`
try:
    subprocess.run(
        [sys.executable, "-m", "panel", "serve", "main.py", "--show"],
        env=env,
        check=True,
    )
except KeyboardInterrupt:
    sys.exit(0)
except subprocess.CalledProcessError as e:
    print(f"ClusterF launch failed: {e}")
    sys.exit(1)
