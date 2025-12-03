import os
import sys
from pathlib import Path

def create_project_dump(source_directory, output_filename, included_extensions=None, ignore_dirs=None):
    # 1. Define standard garbage directories we almost always want to ignore
    files_to_ignore = {'.git', '__pycache__', '.idea', '.vscode', 'node_modules', 
                       'venv', 'env', '.DS_Store', 'dist', 'build', 'coverage', 'target'}
    
    # 2. If user provided specific directories, add them to the set
    if ignore_dirs:
        files_to_ignore.update(ignore_dirs)
    
    always_ignore_extensions = {'.png', '.jpg', '.jpeg', '.gif', '.pyc', '.exe', '.bin', 
                                '.svg', '.ico', '.zip', '.pdf', '.woff', '.woff2', '.ttf', 
                                '.eot', '.mp4', '.mp3', '.wav'}

    files_to_process = []

    # --- Step 1: Walk the directory ---
    for root, dirs, files in os.walk(source_directory):
        # ---------------------------------------------------------
        # This line is the magic: Modifying 'dirs' in-place tells 
        # os.walk to NOT traverse into these directories.
        # This excludes the directory AND all files inside it.
        # ---------------------------------------------------------
        dirs[:] = [d for d in dirs if d not in files_to_ignore]
        
        for file in files:
            file_extension = os.path.splitext(file)[1].lower()
            
            if file_extension in always_ignore_extensions:
                continue

            if included_extensions:
                if file_extension not in included_extensions:
                    continue
                
            full_path = os.path.join(root, file)
            relative_path = os.path.relpath(full_path, source_directory)
            
            if os.path.abspath(output_filename) == os.path.abspath(full_path):
                continue

            files_to_process.append((full_path, relative_path))

    files_to_process.sort(key=lambda x: x[1])

    try:
        with open(output_filename, 'w', encoding='utf-8') as outfile:
            outfile.write(f"\n\n\n\n\n\n\n\n\n\n\n# Project Dump\n")
            outfile.write("\n\n# Directory Structure:\n\n")
            
            if not files_to_process:
                outfile.write("No files found matching criteria.\n")
            
            for _, relative_path in files_to_process:
                outfile.write(f"- {relative_path}\n")
            
            outfile.write("\n\n# Files Content\n\n")

            for full_path, relative_path in files_to_process:
                outfile.write(f"--- START OF FILE: {relative_path} ---\n")
                outfile.write("```\n")
                
                try:
                    with open(full_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                        outfile.write(content)
                except UnicodeDecodeError:
                    outfile.write("[ERROR: Binary or non-UTF-8 file.]")
                except Exception as e:
                    outfile.write(f"[ERROR: {e}]")
                
                outfile.write("\n```\n")
                outfile.write(f"--- END OF FILE: {relative_path} ---\n\n")

            outfile.write(f"# Ignored Directories: {', '.join(sorted(list(files_to_ignore)))}\n")
                
        print(f"Success! Processed {len(files_to_process)} files.")
        print(f"File saved to: {os.path.abspath(output_filename)}")

    except Exception as e:
        print(f"An error occurred: {e}")

def resolve_extensions(ext_input):
    """Helper to convert string inputs like 'web' into lists"""
    if not ext_input:
        return None
    
    if isinstance(ext_input, list):
        return ext_input

    web_ext = ['.md', '.html', '.htm', '.css', '.scss', '.js', '.jsx', '.ts', '.tsx', '.php', '.json', '.vue']
    python_ext = ['.py', '.md']

    if ext_input == 'web':
        return web_ext
    elif ext_input == 'web_py':
        return web_ext + python_ext
    
    # Handle string input like ".py .js"
    ext_list = []
    for ext in ext_input.split():
        if not ext.startswith('.'):
            ext = '.' + ext
        ext_list.append(ext.lower())
    return ext_list

# --- Configuration & Execution ---
if __name__ == "__main__":
    
    # ==========================================
    # CONFIGURATION
    # ==========================================
    
    # 1. Source Directory
    SOURCE_DIR = str(next((p for p in Path(__file__).resolve().parents if (p / "assets").exists()), Path.cwd()))

    # 2. Extensions
    PRESET_EXTENSIONS = 'web_py' 

    # 3. Output Filename
    OUTPUT_FILE = "project_context.txt"

    # 4. CUSTOM IGNORE DIRECTORIES (NEW FEATURE)
    #    Add the folder names you want to completely skip here.
    #    Example: {'tests', 'docs', 'legacy_code'}
    CUSTOM_IGNORES = {'data', 'tmp', 'logs', 'exports'} 

    # ==========================================
    # LOGIC
    # ==========================================

    final_dir = SOURCE_DIR
    final_exts = None

    # --- Mode A: Hardcoded (Silent Run) ---
    if SOURCE_DIR:
        if not os.path.exists(SOURCE_DIR):
            print(f"Error: Hardcoded directory not found: {SOURCE_DIR}")
            sys.exit(1)
        
        final_exts = resolve_extensions(PRESET_EXTENSIONS)
        print(f"--- Processing Hardcoded Path: {SOURCE_DIR} ---")
        print(f"--- Ignoring Directories: {CUSTOM_IGNORES} ---")

    # --- Mode B: Interactive (Ask User) ---
    else:
        print("--- Directory to Text Converter ---")
        
        # 1. Ask for Directory
        while True:
            raw_path = input("Enter the full path (or press ENTER for Parent Directory):\n> ").strip()
            
            if not raw_path:
                parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
                print(f"\nYou selected the PARENT directory: {parent_dir}")
                if input("Are you sure? (y/n): ").strip().lower() == 'y':
                    final_dir = parent_dir
                    break
                continue

            if raw_path.startswith(('"', "'")) and raw_path.endswith(('"', "'")):
                raw_path = raw_path[1:-1]
                
            if os.path.isdir(raw_path):
                final_dir = raw_path
                break
            else:
                print(f"Error: Directory not found at '{raw_path}'.\n")

        # 2. Ask for Extensions
        print("\nWhich extensions? (ENTER for all, 'web' for web files, or type '.py .js')")
        ext_input = input("> ").strip().lower()
        final_exts = resolve_extensions(ext_input)
        
        print("\n--- Processing ---")

    # --- Run with Custom Ignores ---
    create_project_dump(
        final_dir, 
        OUTPUT_FILE, 
        included_extensions=final_exts, 
        ignore_dirs=CUSTOM_IGNORES  # <--- Passing your specific ignores here
    )
