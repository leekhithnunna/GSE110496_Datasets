import os
import shutil

def create_processed_folder(base_path: str) -> str:
    """Creates a 'processed' folder in the given base path if it doesn't exist."""
    processed_folder = os.path.join(base_path, "processed")
    os.makedirs(processed_folder, exist_ok=True)
    return processed_folder

def rename_folder_without_extension(folder_path: str, extension: str) -> str:
    """Renames a folder by removing the specified extension."""
    new_folder_name = folder_path[:-len(extension)]
    os.rename(folder_path, new_folder_name)
    return new_folder_name

def process_folder(folder_path: str, extension: str = ".tsv"):
    """Processes folders by removing the extension and moving them to the 'processed' directory."""
    if not os.path.exists(folder_path):
        print(f"The folder '{folder_path}' does not exist.")
        return
    
    processed_folder = create_processed_folder(folder_path)

    for item_name in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item_name)

        if os.path.isdir(item_path) and item_name.endswith(extension):
            try:
                new_folder_name = rename_folder_without_extension(item_path, extension)
                new_folder_path = os.path.join(folder_path, new_folder_name)
                
                shutil.move(new_folder_path, processed_folder)
                print(f"Renamed and moved: {item_name} -> {new_folder_name}")
            except Exception as e:
                print(f"Error processing {item_name}: {e}")
        else:
            print(f"Skipped: {item_name}")

# Example usage
if __name__ == "__main__":
    folder_path = ""
    process_folder(folder_path)
