import os


if __name__ == "__main__":
    folder = "../experiment/graphs"
    for filename in os.listdir(folder):
        path = os.path.join(folder, filename)
        print("./local " + path)
