from pathlib import Path


SCHEMA = "date-{date}_prot-{prot}_driver-{driver}_cons-{cons}_exp-{exp}_emb-{i}"
FOLDER = Path("data", "enrichment")


if __name__ == "__main__":

    i = 5

    for file in FOLDER.iterdir():
        if not file.is_file() or file.name.startswith("_"):
            continue

        driver, cons, exp = "endo", "pdzb3", "homo"
        date, _, prot, *_ = file.stem.split("_")
        date = date.replace("-", "")

        name = SCHEMA.format(
            date=date, prot=prot, driver=driver, cons=cons, exp=exp, i=i
        )

        if not (path := FOLDER.joinpath(name)).exists():
            path.mkdir()
            file.rename(path.joinpath(name).with_suffix(".tif"))

        i += 1
