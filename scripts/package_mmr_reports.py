#!/usr/bin/env python3
import os
import shutil
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Package clinical MMR reports, plots, and dashboards into a single shareable ZIP file.")
    parser.add_argument("--results-dir", required=True, help="Base results directory containing sample folders.")
    parser.add_argument("--output-zip", required=True, help="Path to write the output ZIP file (e.g. cohort_reports.zip).")
    return parser.parse_args()

def main():
    args = parse_args()
    results_dir = os.path.abspath(args.results_dir)
    output_zip = os.path.abspath(args.output_zip)
    
    if not os.path.exists(results_dir):
        print(f"Error: Results directory does not exist: {results_dir}")
        return
        
    # Get scripts directory relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    compile_script = os.path.join(script_dir, "compile_mmr_summary.py")
    
    # 1. Compile cohort dashboard to make sure it's up to date
    print("Compiling cohort dashboard...")
    if os.path.exists(compile_script):
        cmd = f"python3 {compile_script} --results-dir {results_dir}"
        os.system(cmd)
    else:
        print("Warning: compile_mmr_summary.py not found in same directory. Skipping regeneration.")
        
    dashboard_path = os.path.join(results_dir, "cohort_mmr_dashboard.html")
    if not os.path.exists(dashboard_path):
        print("Error: cohort_mmr_dashboard.html not found. Please run compilation first.")
        return
        
    # 2. Create temporary packaging directory
    temp_dir = os.path.join(results_dir, "mmr_reports_package_temp")
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir, ignore_errors=True)
    os.makedirs(temp_dir)
    
    # 3. Copy dashboard
    shutil.copy2(dashboard_path, os.path.join(temp_dir, "cohort_mmr_dashboard.html"))
    
    # 4. Copy summary subdirectories
    print("Extracting clinical report assets (HTMLs, BAF plots, IGV reports)...")
    tar_cmd = f"cd {results_dir} && tar -cf - */summary | (cd {temp_dir} && tar -xf -)"
    os.system(tar_cmd)
    
    # 5. Create ZIP archive
    print(f"Creating ZIP archive: {output_zip} ...")
    os.system(f"cd {temp_dir} && zip -q -r {output_zip} .")
    
    # 6. Clean up temporary directory
    shutil.rmtree(temp_dir, ignore_errors=True)
    print("Packaging complete! You can download and share the resulting ZIP file.")

if __name__ == "__main__":
    main()
