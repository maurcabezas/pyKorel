import os
import sys
import subprocess
import configparser
import zipfile
import shutil

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
    extract_dir = "etc/spec_test"

    # Check if spec_test.zip exists
    if not os.path.isfile(zip_path):
        print(f"Warning: {zip_path} not found in the repository root. Skipping extraction.")
        print("Please ensure spec_test.zip is present or manually create the spec_test directory with test FITS files.")
        return

    try:
        # Create spec_test directory if it doesn't exist
        os.makedirs(extract_dir, exist_ok=True)
        
        # Extract the zip file
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        print(f"Successfully extracted {zip_path} to {extract_dir}/")

        # Verify extraction
        if not os.listdir(extract_dir):
            print(f"Error: {extract_dir} is empty after extraction. The zip file may be corrupted.")
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
        if package.strip() and not package.startswith('#'):  # Skip empty lines and comments
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
    """Update initpath in relevant Python files (placeholder for future use)."""
    # Currently, only main.cfg is updated. Add logic here if other files need initpath updates.
    print("No additional files require initpath updates at this time.")

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
    print("  cd spec_test")
    print("  python ../src/prekor.py")
