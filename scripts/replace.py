import os
from pathlib import Path
from pathlib import PurePath
import shutil


def allfiles(path: Path):
    for entry in path.iterdir():
        if entry.is_file():
            yield entry
        elif entry.is_dir():
            yield from allfiles(entry)
        else:
            yield None


if __name__ == "__main__":

    OLD: str = "./fortune_business_theme_from_old_server"
    NEW: str = "./fortune_business_theme_d9"

    for folder in Path(OLD).iterdir():
        if folder.name not in ["css", "templates"]:
            continue

        for file in allfiles(folder):
            src = file.path
            dst = Path(NEW, *src.parts[1:])

            print(src)
            print(dst)
