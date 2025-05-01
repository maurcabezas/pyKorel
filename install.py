import os
import sys
import subprocess
import configparser
import zipfile
import re

def get_current_path():
    """Return the current working directory."""
    return os.getcwd()

def check_python_version():
    """Ensure Python version is 3.8 or higher."""
    required_version = (3, 8)
    current_version = sys.version_info[:2]
    if current_version < required_version:
        print(f"Error: Python {required_version[0]}.{required_version[1]} or higher is required. Found {current_version[0]}.{current_version[1]}.")
        sys.exit(1)
    print(f"Python version {current_version[0]}.{current_version[1]} is compatible.")

def unzip_spec_test():
    """Unzip spec_test.zip to the spec_test directory."""
    zip_path = "etc/spec_test.zip"
    extract_dir = "etc/"

    if not os.path.isfile(zip_path):
        print(f"Warning: {zip_path} not found in the repository root. Skipping extraction.")
        print("Please ensure spec_test.zip is present or manually create the spec_test directory with test FITS files.")
        return

    try:
        os.makedirs(extract_dir, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        print(f"Successfully extracted {zip_path} to {extract_dir}/")
        if not os.listdir(extract_dir):
            print(f"Warning: {extract_dir} is empty after extraction. The zip file may be corrupted.")
        else:
            print(f"Contents of {extract_dir}: {os.listdir(extract_dir)}")
    except zipfile.BadZipFile:
        print(f"Error: {zip_path} is corrupted or not a valid zip file.")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when extracting to {extract_dir}. Run the script with appropriate permissions.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to extract {zip_path}. Details: {e}")
        sys.exit(1)

def check_and_install_requirements():
    """Install dependencies listed in requirements.txt."""
    requirements_file = "requirements.txt"
    
    if not os.path.isfile(requirements_file):
        print(f"Error: {requirements_file} not found in the repository root.")
        sys.exit(1)

    with open(requirements_file, 'r') as f:
        requirements = f.read().splitlines()

    print("Installing dependencies...")
    for package in requirements:
        if package.strip() and not package.startswith('#'):
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", package])
                print(f"Successfully installed {package}")
            except subprocess.CalledProcessError as e:
                print(f"Failed to install {package}. Error: {e}")
                sys.exit(1)

def update_main_cfg(path):
    """Update the initpath in src/main.cfg with the current path."""
    config_file = "src/main.cfg"
    
    if not os.path.isfile(config_file):
        print(f"Error: {config_file} not found.")
        sys.exit(1)

    try:
        config = configparser.ConfigParser()
        config.read(config_file)
        if 'pyKorel config' not in config:
            config['pyKorel config'] = {}
        config['pyKorel config']['initpath'] = path

        with open(config_file, 'w') as configfile:
            config.write(configfile)
        print(f"Updated {config_file} with initpath: {path}")
    except Exception as e:
        print(f"Error updating {config_file}: {e}")
        sys.exit(1)

def update_initpath_in_files(path):
    """Update initpath in specified Python files."""
    files = ['src/prekor.py', 'src/multiprekor.py']  # Add 'src/korel.py', 'src/multikorel.py' if needed
    pattern = re.compile(r"initpath\s*=\s*['\"].*?['\"]")

    for file in files:
        if not os.path.isfile(file):
            print(f"Warning: {file} not found. Skipping initpath update for this file.")
            continue

        try:
            with open(file, 'r') as f:
                content = f.read()

            # Replace the initpath line
            new_content = pattern.sub(f"initpath = '{path}'", content)
            if content == new_content:
                print(f"Warning: No initpath line found in {file}. Ensure it contains 'initpath = ...'.")
                continue

            with open(file, 'w') as f:
                f.write(new_content)
            print(f"Updated initpath in {file} to: {path}")
        except PermissionError:
            print(f"Error: Permission denied when updating {file}. Run the script with appropriate permissions.")
            sys.exit(1)
        except Exception as e:
            print(f"Error updating {file}: {e}")
            sys.exit(1)

if __name__ == "__main__":
    print("Starting pyKorel installation...")

    # Check Python version
    check_python_version()

    # Unzip spec_test.zip
    unzip_spec_test()

    # Get current path
    current_path = get_current_path()
    print(f"Current working directory: {current_path}")

    # Update configuration and files
    update_initpath_in_files(current_path)
    update_main_cfg(current_path)

    # Install dependencies
    check_and_install_requirements()

    print("pyKorel installation completed successfully!")
    print("To test the installation, navigate to the spec_test directory and run:")
    print("  cd etc/spec_test")
    print("  python ../../src/prekor.py")
