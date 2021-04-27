import gzip
import os
import requests

LONGJING43 = "longjing_43_transcriptome.fas"
TEA015198 = "tea015198.fas"

DATA_URLS = {
    LONGJING43: {
        "url": "http://tpia.teaplant.org/web/Download/Transcriptome_data/Camellia_sinensis_var_sinensis_cv_Longjing_43_Trans.fas.gz",
        "compression": "gzip",
    },
    TEA015198: {
        "url": "http://tpia.teaplant.org/Download_Sequence?ID=TEA015198",
    },
}


def download_data():
    if not os.path.exists("./data"):
        os.mkdir("./data")
    for key in DATA_URLS.keys():
        res = requests.get(DATA_URLS[key]["url"])
        data = res.content
        if DATA_URLS[key].get("compression") == "gzip":
            data = gzip.decompress(data)
        with open("./data/" + key, "wb") as f:
            f.write(data)


if __name__ == "__main__":
    download_data()
