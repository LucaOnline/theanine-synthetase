import os

import requests

TEA015198 = "tea015198.fasta"

DATA_URLS = {TEA015198: "http://tpia.teaplant.org/Download_Sequence?ID=TEA015198"}


def download_data():
    res = requests.get(DATA_URLS[TEA015198])
    if not os.path.exists("./data"):
        os.mkdir("./data")
    with open("./data/" + TEA015198, "wb") as f:
        f.write(bytes(res.text, res.encoding))
