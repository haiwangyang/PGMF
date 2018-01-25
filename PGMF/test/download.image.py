import random
import urllib.request

def download_web_image(url):
    name = "../data/output/" + str(random.randrange(1,1000)) + ".jpg"
    urllib.request.urlretrieve(url, name)

if __name__ == '__main__': 
    download_web_image("https://is1-ssl.mzstatic.com/image/thumb/Purple118/v4/19/4a/21/194a21e8-2eeb-46a8-2ec9-632aae23b1a3/AppIcon-1x_U007emarketing-85-220-0-5.png/200x0w.jpg")

    
