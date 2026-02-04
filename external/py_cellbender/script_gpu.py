import os
import subprocess

# --- Hardware Configuration ---
USE_GPU = True  # Toggle: True for --cuda, False for CPU
 
with open("input.txt", "r") as file:
    INPUT_H5 = file.readline().strip()

if not INPUT_H5:
    raise ValueError("input.txt is empty")

if not os.path.exists(INPUT_H5):
    raise FileNotFoundError(f"Input file not found: {INPUT_H5}")
    
print("Input file (raw_feature_bc_matrix.h5):", INPUT_H5)

output_h5 = "output.h5"
temp_checkpoint = "ckpt.tar.gz"
final_checkpoint = "final_checkpoint.tar.gz"

cmd = [
    "cellbender", "remove-background", 
    "--input", INPUT_H5,
    "--output", output_h5,
    "--expected-cells", "10000",
    "--epochs", "150", 
    "--fpr", "0.01",
    "--learning-rate", "5e-5",
    "--total-droplets-included", "40000"
]

if USE_GPU:
    cmd.append("--cuda")

# --- Execution ---
print(f"Starting CellBender")
if not os.path.exists(INPUT_H5):
    raise FileNotFoundError(f"Input file not found: {INPUT_H5}")

try:
    result = subprocess.run(
        cmd, 
        check=True, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True        
    )
    
    # Rename checkpoint if it was created
    #if os.path.exists(temp_checkpoint):
    #    os.rename(temp_checkpoint, final_checkpoint)
    #    print(f"Checkpoint saved to: {final_checkpoint}")

    print(f"CellBender complete. Output: {output_h5}")

except subprocess.CalledProcessError as e:
    print(f"CellBender failed.")
    print(f"STDOUT: {e.stdout}")
    print(f"STDERR: {e.stderr}")