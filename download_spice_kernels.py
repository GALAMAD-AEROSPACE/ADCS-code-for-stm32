import os
import requests
from tqdm import tqdm
import sys

def download_file(url, filename):
    """
    Download a file with progress bar
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        # Get total file size
        total_size = int(response.headers.get('content-length', 0))
        
        # Create progress bar
        progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True, desc=filename)
        
        # Download the file
        with open(filename, 'wb') as f:
            for data in response.iter_content(chunk_size=8192):
                progress_bar.update(len(data))
                f.write(data)
        
        progress_bar.close()
        return True
    except Exception as e:
        print(f"Error downloading {filename}: {str(e)}")
        return False

def main():
    # Create kernels directory if it doesn't exist
    if not os.path.exists('kernels'):
        os.makedirs('kernels')
    
    # Base URL for NAIF kernels
    base_url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels"
    
    # Define the kernels to download
    kernels = {
        'leapseconds.tls': f"{base_url}/lsk/latest_leapseconds.tls",
        'pck00010.tpc': f"{base_url}/pck/pck00010.tpc",
        'de430.bsp': f"{base_url}/spk/planets/de430.bsp"  # Using de430 instead of de440 for smaller size
    }
    
    print("Starting download of SPICE kernels...")
    
    # Download each kernel
    for filename, url in kernels.items():
        print(f"\nDownloading {filename}...")
        if download_file(url, os.path.join('kernels', filename)):
            print(f"Successfully downloaded {filename}")
        else:
            print(f"Failed to download {filename}")
            sys.exit(1)
    
    print("\nAll kernels downloaded successfully!")
    print("Files are located in the 'kernels' directory")
    print("\nNote: You'll need to update the paths in your C++ code to point to these files.")
    print("For example, change 'leapseconds.tls' to 'kernels/leapseconds.tls'")

if __name__ == "__main__":
    main() 