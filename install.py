import os
import re
import configparser
import subprocess
import sys

def update_initpath_in_files(path):
    files = ['src/prekor.py', 'src/multiprekor.py']#, 'src/korel.py', 'src/multikorel.py']
    
    # Define the regex pattern to match the initpath line
    pattern = re.compile(r"initpath\s*=\s*['\"].*?['\"]")
    # pattern = re.match(r'^initpath', line)

    for file in files:
        with open(file, 'r') as f:
            content = f.read()
        
        # Replace the initpath line with the new path
        new_content = pattern.sub(f"initpath = '{path}'", content)
        
        with open(file, 'w') as f:
            f.write(new_content)
        
        print(f"Updated initpath in {file}")


def get_current_path():
    return os.getcwd()

def check_and_install_requirements():
    with open('requirements.txt', 'r') as f:
        requirements = f.read().splitlines()
    
    for package in requirements:
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        except subprocess.CalledProcessError as e:
            print(f"Failed to install {package}. Error: {e}")

def update_main_cfg(path):
    config = configparser.ConfigParser()
    config.read('src/main.cfg')
    config['pyKorel config']['initpath'] = path
    
    with open('src/main.cfg', 'w') as configfile:
        config.write(configfile)
    print("Updated main.cfg with the new path")

if __name__ == "__main__":
    current_path = get_current_path()
    update_initpath_in_files(current_path)
    update_main_cfg(current_path)
    check_and_install_requirements()
