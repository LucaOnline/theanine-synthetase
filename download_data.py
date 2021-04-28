import gzip
import os
import requests

from data_index import DATA_URLS


def download_data():
    if not os.path.exists("./data"):
        os.mkdir("./data")
    for key in DATA_URLS.keys():
        if os.path.isfile("./data/" + key):
            continue
        res = requests.get(DATA_URLS[key]["url"])
        data = res.content
        if DATA_URLS[key].get("compression") == "gzip":
            data = gzip.decompress(data)
        with open("./data/" + key, "wb") as f:
            f.write(data)


if __name__ == "__main__":
    download_data()
