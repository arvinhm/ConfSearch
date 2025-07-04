#!/usr/bin/env python3
"""
Automated MacroModel Conformational Search Runner
Automatically runs MacroModel with your .mae and .com files
"""

import os
import subprocess
import sys
import time
from pathlib import Path

class MacroModelRunner:
    """Automated MacroModel conformational search runner"""
    
    def __init__(self):
        self.schrodinger_path = os.environ.get('SCHRODINGER', '')
        if not self.schrodinger_path:
            # Try common Windows installation paths
            possible_paths = [
                r"C:\Program Files\Schrodinger2025-2",
                r"C:\Program Files\Schrodinger2024-4",
                r"C:\Program Files\Schrodinger2024-3",
                r"C:\Program Files\Schrodinger2024-2",
                r"C:\Program Files\Schrodinger2024-1"
            ]
            
            for path in possible_paths:
                if os.path.exists(path):
                    self.schrodinger_path = path
                    break
        
        if not self.schrodinger_path:
            raise RuntimeError("Schrödinger installation not found! Please set SCHRODINGER environment variable.")
        
        self.macromodel_exe = os.path.join(self.schrodinger_path, "macromodel.exe")
        if not os.path.exists(self.macromodel_exe):
            self.macromodel_exe = "macromodel"  # Try from PATH

    def run_conformational_search(self, com_file, wait=True, verbose=True):
        """
        Run MacroModel conformational search
        
        Args:
            com_file: Path to .com file
            wait: Wait for job to complete
            verbose: Print progress information
        """
        
        com_path = Path(com_file)
        if not com_path.exists():
            raise FileNotFoundError(f"COM file not found: {com_file}")
        
        # Check if corresponding .mae file exists
        mae_files = list(com_path.parent.glob("*.mae"))
        if not mae_files:
            print("⚠️  Warning: No .mae files found in directory")
        
        if verbose:
            print(f"🚀 Starting MacroModel conformational search...")
            print(f"📁 Working directory: {com_path.parent}")
            print(f"📄 COM file: {com_path.name}")
            print(f"🔧 Schrödinger path: {self.schrodinger_path}")
        
        # Build command
        if wait:
            cmd = [self.macromodel_exe, "-WAIT", str(com_path)]
        else:
            cmd = [self.macromodel_exe, str(com_path)]
        
        if verbose:
            print(f"⚡ Command: {' '.join(cmd)}")
            print("=" * 60)
        
        try:
            # Change to the directory containing the files
            original_dir = os.getcwd()
            os.chdir(com_path.parent)
            
            # Run MacroModel
            if verbose:
                # Run with real-time output
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                    bufsize=1
                )
                
                print("📊 MacroModel Output:")
                print("-" * 40)
                
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        print(output.strip())
                
                return_code = process.poll()
                
            else:
                # Run silently
                result = subprocess.run(cmd, capture_output=True, text=True)
                return_code = result.returncode
                if return_code != 0:
                    print(f"❌ Error: {result.stderr}")
            
            # Change back to original directory
            os.chdir(original_dir)
            
            if return_code == 0:
                if verbose:
                    print("\n" + "=" * 60)
                    print("🎉 MacroModel job completed successfully!")
                    self._show_output_files(com_path)
                return True
            else:
                print(f"❌ MacroModel job failed with return code: {return_code}")
                return False
                
        except FileNotFoundError:
            print(f"❌ MacroModel executable not found: {self.macromodel_exe}")
            print("💡 Make sure Schrödinger is installed and SCHRODINGER environment variable is set")
            return False
        except Exception as e:
            print(f"❌ Error running MacroModel: {e}")
            return False

    def _show_output_files(self, com_path):
        """Show the output files created"""
        base_name = com_path.stem
        output_dir = com_path.parent
        
        # Look for output files
        output_extensions = ['.maegz', '.mae', '.log', '.out', '.tmp']
        output_files = []
        
        for ext in output_extensions:
            pattern = f"*{ext}"
            files = list(output_dir.glob(pattern))
            output_files.extend(files)
        
        if output_files:
            print("\n📁 Output files created:")
            for file in sorted(output_files):
                size_mb = file.stat().st_size / (1024 * 1024)
                print(f"   ✅ {file.name} ({size_mb:.2f} MB)")
        else:
            print("\n⚠️  No output files found")

    def batch_run(self, directory=".", pattern="*.com", max_concurrent=2):
        """
        Run multiple conformational searches in batch
        
        Args:
            directory: Directory containing .com files
            pattern: Pattern to match .com files
            max_concurrent: Maximum number of concurrent jobs
        """
        
        dir_path = Path(directory)
        com_files = list(dir_path.glob(pattern))
        
        if not com_files:
            print(f"❌ No .com files found matching pattern: {pattern}")
            return
        
        print(f"🎯 Found {len(com_files)} .com files to process")
        print(f"🔄 Running {max_concurrent} jobs concurrently")
        print("=" * 60)
        
        # Run jobs in batches
        running_jobs = []
        completed = 0
        
        for i, com_file in enumerate(com_files):
            print(f"\n📋 Job {i+1}/{len(com_files)}: {com_file.name}")
            
            # Start job without waiting
            cmd = [self.macromodel_exe, str(com_file)]
            
            try:
                original_dir = os.getcwd()
                os.chdir(com_file.parent)
                
                process = subprocess.Popen(cmd)
                running_jobs.append((process, com_file.name))
                
                os.chdir(original_dir)
                
                print(f"✅ Started job: {com_file.name}")
                
                # Wait if we've reached max concurrent jobs
                if len(running_jobs) >= max_concurrent:
                    self._wait_for_job_completion(running_jobs)
                    completed += 1
                    
            except Exception as e:
                print(f"❌ Failed to start job {com_file.name}: {e}")
        
        # Wait for remaining jobs
        while running_jobs:
            self._wait_for_job_completion(running_jobs)
            completed += 1
        
        print(f"\n🎉 All {len(com_files)} conformational searches completed!")

    def _wait_for_job_completion(self, running_jobs):
        """Wait for at least one job to complete"""
        while running_jobs:
            for i, (process, name) in enumerate(running_jobs):
                if process.poll() is not None:
                    print(f"✅ Completed: {name}")
                    running_jobs.pop(i)
                    return
            time.sleep(5)  # Check every 5 seconds

def main():
    """Main function with command line interface"""
    
    runner = MacroModelRunner()
    
    print("🧬 MacroModel Conformational Search Runner")
    print("=" * 50)
    
    # Check if specific .com file provided as argument
    if len(sys.argv) > 1:
        com_file = sys.argv[1]
        print(f"Running single job: {com_file}")
        success = runner.run_conformational_search(com_file)
        if success:
            print("🎉 Job completed successfully!")
        else:
            print("❌ Job failed!")
        return
    
    # Interactive mode
    current_dir = Path(".")
    com_files = list(current_dir.glob("*.com"))
    
    if not com_files:
        print("❌ No .com files found in current directory")
        print("💡 Make sure you're in the directory with your .com files")
        return
    
    print(f"📁 Found {len(com_files)} .com files:")
    for i, com_file in enumerate(com_files, 1):
        print(f"   {i}. {com_file.name}")
    
    print("\nOptions:")
    print("1. Run single conformational search")
    print("2. Run all conformational searches (batch)")
    print("3. Exit")
    
    try:
        choice = input("\nEnter your choice (1-3): ").strip()
        
        if choice == "1":
            # Single job
            print("\nSelect .com file:")
            for i, com_file in enumerate(com_files, 1):
                print(f"   {i}. {com_file.name}")
            
            file_choice = int(input("Enter file number: ")) - 1
            if 0 <= file_choice < len(com_files):
                selected_com = com_files[file_choice]
                print(f"\n🚀 Running: {selected_com.name}")
                success = runner.run_conformational_search(selected_com)
                if success:
                    print("🎉 Job completed successfully!")
                else:
                    print("❌ Job failed!")
            else:
                print("❌ Invalid file number")
                
        elif choice == "2":
            # Batch jobs
            max_concurrent = input("Enter max concurrent jobs (default: 2): ").strip()
            if not max_concurrent:
                max_concurrent = 2
            else:
                max_concurrent = int(max_concurrent)
            
            print(f"\n🔄 Running all {len(com_files)} jobs...")
            runner.batch_run(max_concurrent=max_concurrent)
            
        elif choice == "3":
            print("👋 Goodbye!")
            
        else:
            print("❌ Invalid choice")
            
    except KeyboardInterrupt:
        print("\n\n👋 Interrupted by user")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()