TOOLS = [f"{CACHE}/bin/unmux", f"{CACHE}/bin/lima"]


def download(key, name):
    return f"{CACHE}/downloads/{key}/{name}"


def downloads(inputs):
    return [download(key, name) for key, name, _ in inputs]


def work(recipe, name):
    return f"{CACHE}/work/{recipe}/{name}"


def result(name):
    return f"{RESULTS}/{name}"


def logfile(name):
    return f"{RESULTS}/logs/{name}.log"
